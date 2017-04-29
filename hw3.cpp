/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: Aditya Aggarwal
 * *************************
 */

#ifdef WIN32
#include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
#include <GL/gl.h>
#include <GL/glut.h>
#elif defined(__APPLE__)
#include <OpenGL/gl.h>
#include <GLUT/glut.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef WIN32
#define strcasecmp _stricmp
#endif

#include <imageIO.h>
#include <cmath>
#include <tgmath.h>
#include <glm/glm.hpp>
#include <iostream>
#include "glm/ext.hpp"
#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

char * filename = NULL;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 60.0

unsigned char buffer[HEIGHT][WIDTH][3];

struct Vertex
{
    double position[3];
    double color_diffuse[3];
    double color_specular[3];
    double normal[3];
    double shininess;
};

struct Triangle
{
    Vertex v[3];
};

struct Sphere
{
    double position[3];
    double color_diffuse[3];
    double color_specular[3];
    double shininess;
    double radius;
};

struct Light
{
    double position[3];
    double color[3];
};

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);


//triangle inside outside test
glm::vec3 barycentric_test(float parameter, Triangle triangle, glm::vec3 direction, glm::vec3 initial) {
    
    //location of potential intersection
    glm::vec3 location = (direction * parameter) + initial;
    
    
    
    //vertex positions
    glm::vec3 vertexa = glm::vec3(triangle.v[0].position[0], triangle.v[0].position[1], triangle.v[0].position[2]);
    glm::vec3 vertexb = glm::vec3(triangle.v[1].position[0], triangle.v[1].position[1], triangle.v[1].position[2]);
    glm::vec3 vertexc = glm::vec3(triangle.v[2].position[0], triangle.v[2].position[1], triangle.v[2].position[2]);
    
    
    
    //project into 2D
    //ignore z, x, or y dimension?
    bool xPlane = !((vertexa.x != vertexb.x) || (vertexb.x != vertexc.x));
    bool yPlane = !((vertexa.y != vertexb.y) || (vertexb.y != vertexc.y));
    bool zPlane = !((vertexa.z != vertexb.z) || (vertexb.z != vertexc.z));
    
    //declaration of barycentric coordinates
    float alpha;
    float beta;
    float gamma;
    
    if(!xPlane && !yPlane) { //xy plane projection
        
        float mainTriangleArea = ((vertexb.x - vertexa.x)*(vertexc.y- vertexa.y) - (vertexc.x - vertexa.x)*(vertexb.y-vertexa.y));
        mainTriangleArea *= 0.5;
        alpha = ((vertexb.x - location.x)*(vertexc.y-location.y) - (vertexc.x - location.x)*(vertexb.y-location.y));
        alpha *= 0.5;
        beta = ((location.x - vertexa.x)*(vertexc.y-vertexa.y) - (vertexc.x-vertexa.x)*(location.y - vertexa.y));
        beta *= 0.5;
        gamma = ((vertexb.x - vertexa.x)*(location.y-vertexa.y) - (location.x-vertexa.x)*(vertexb.y - vertexa.y));
        gamma *= 0.5;
        
        
        alpha/=mainTriangleArea;
        beta/=mainTriangleArea;
        gamma/=mainTriangleArea;
        
        
        
    } else if (!xPlane && !zPlane) { //xz plane projection
        float mainTriangleArea = ((vertexb.x - vertexa.x)*(vertexc.z- vertexa.z) - (vertexc.x - vertexa.x)*(vertexb.z-vertexa.z));
        mainTriangleArea *=0.5;
        alpha = ((vertexb.x - location.x)*(vertexc.z-location.z) - (vertexc.x - location.x)*(vertexb.z-location.z));
        alpha *= 0.5;
        beta = ((location.x - vertexa.x)*(vertexc.z-vertexa.z) - (vertexc.x-vertexa.x)*(location.z - vertexa.z));
        beta *= 0.5;
        gamma = ((vertexb.x - vertexa.x)*(location.z-vertexa.z) - (location.x-vertexa.x)*(vertexb.z - vertexa.z));
        gamma *= 0.5;
        
        
        alpha/=mainTriangleArea;
        beta/=mainTriangleArea;
        gamma/=mainTriangleArea;
        
    } else { //yz plane projection
        float mainTriangleArea = ((vertexb.y - vertexa.y)*(vertexc.z- vertexa.z) - (vertexc.y - vertexa.y)*(vertexb.z-vertexa.z));
        mainTriangleArea *=0.5;
        alpha = ((vertexb.y - location.y)*(vertexc.z-location.z) - (vertexc.y - location.y)*(vertexb.z-location.z));
        alpha *= 0.5;
        beta = ((location.y - vertexa.y)*(vertexc.z-vertexa.z) - (vertexc.y-vertexa.y)*(location.z - vertexa.z));
        beta *= 0.5;
        gamma = ((vertexb.y - vertexa.y)*(location.z-vertexa.z) - (location.y-vertexa.y)*(vertexb.z - vertexa.z));
        gamma *= 0.5;
        
        
        alpha/=mainTriangleArea;
        beta/=mainTriangleArea;
        gamma/=mainTriangleArea;
        
        
    }
    
    
    return glm::vec3(alpha, beta, gamma);
    
}


std::pair<glm::vec3, int> intersectCalculations(glm::vec3 direction, glm::vec3 initial, Triangle & minTriangle, Sphere & smallestSphere, bool intersectionTest) { //for shadow rays
    //intersection calculations
    float min_intersection;
    float intersectionParam;
    bool hasIntersection = false;
    
    //sphere intersection
    bool hasSphereIntersection = false;
    
    for(int i = 0; i< num_spheres; ++i){
        
        float b = 2 * (direction.x*(initial.x - spheres[i].position[0]) + direction.y*(initial.y - spheres[i].position[1]) + direction.z*(initial.z - spheres[i].position[2]));
        float c = (initial.x - spheres[i].position[0]) * (initial.x - spheres[i].position[0]) + (initial.y - spheres[i].position[1]) * (initial.y - spheres[i].position[1]) + (initial.z - spheres[i].position[2])*(initial.z - spheres[i].position[2]) - (spheres[i].radius * spheres[i].radius);
        float discriminant = (b*b) - 4*c;
        if(discriminant>=0) {
            float discRoot = sqrt(discriminant);
            float roota = ((-1)*b + discRoot)/(2);
            float rootb = ((-1)*b - discRoot)/(2);
            intersectionParam = fmin(roota, rootb);
            
            
           
            
            if(intersectionParam>0) {
                if(!hasIntersection) {
                    hasIntersection = true;
                    hasSphereIntersection = true;
                    min_intersection = intersectionParam;
                    smallestSphere = spheres[i];
                } else if(intersectionParam < min_intersection) {
                    min_intersection = intersectionParam;
                    smallestSphere = spheres[i];
                }
            }
        }
        
    }
    
    
    
    
    
    
    
    
    //triangle intersection
    glm::vec3 barycentric_coordinates = glm::vec3(-1, -1, -1);
    glm::vec3 barycentric_coordinates_temp = glm::vec3(-1, -1, -1);
    bool isTriangleIntersection = false;
    
    
    
    for(int i = 0; i<num_triangles; ++i) {
        glm::vec3 first_vec = glm::vec3(triangles[i].v[0].position[0] - triangles[i].v[1].position[0], triangles[i].v[0].position[1] - triangles[i].v[1].position[1], triangles[i].v[0].position[2] - triangles[i].v[1].position[2]);
        glm::vec3 second_vec = glm::vec3(triangles[i].v[1].position[0] - triangles[i].v[2].position[0], triangles[i].v[1].position[1] - triangles[i].v[2].position[1], triangles[i].v[1].position[2] - triangles[i].v[2].position[2]);
        
        glm::vec3 orthogonal = glm::cross(first_vec, second_vec);
        glm::vec3 normal = glm::normalize(orthogonal);
        glm::vec3 vertex = glm::vec3(triangles[i].v[0].position[0], triangles[i].v[0].position[1], triangles[i].v[0].position[2]);
        
        float plane_d = glm::dot(normal, vertex); //for d
        if(intersectionTest)
            plane_d -= 0.999*glm::dot(normal, initial); // 0.999 makes small value added for shadows proportional to size of magnitude, better than constant to remove artifacts due to floating pt operations
        else
            plane_d -= glm::dot(normal, initial);
        
        float denom = glm::dot(direction, normal);
        if (denom !=0) { //if ray not parallel to plane
            intersectionParam = plane_d/denom;
            
            
            if(intersectionParam>0) {
                if(!hasIntersection){ //if there was no intersection with a sphere
                    barycentric_coordinates_temp = barycentric_test(intersectionParam, triangles[i], direction, initial);
                    if(barycentric_coordinates_temp.x >=0  && barycentric_coordinates_temp.y >= 0 && barycentric_coordinates_temp.z >= 0 && glm::length(barycentric_coordinates_temp) < 1.1) { //for triangle inside outside test
                        isTriangleIntersection = true;
                        hasIntersection = true;
                        min_intersection = intersectionParam;
                        minTriangle = triangles[i];
                        barycentric_coordinates = barycentric_coordinates_temp;
                    }
                } else  if(intersectionParam<min_intersection) { //if there was an intersection with a sphere compare with smallest of them
                    barycentric_coordinates_temp = barycentric_test(intersectionParam, triangles[i], direction, initial);
                    if(barycentric_coordinates_temp.x >=0  && barycentric_coordinates_temp.y >= 0 && barycentric_coordinates_temp.z >= 0 && glm::length(barycentric_coordinates_temp) < 1.1) { //for triangle inside outside test
                        if(!isTriangleIntersection)
                            isTriangleIntersection = true;
                        min_intersection = intersectionParam;
                        minTriangle = triangles[i];
                        barycentric_coordinates = barycentric_coordinates_temp;
                    }
                }
            }
        }
        
    }
    
    
    
    if(hasIntersection) {
        if(isTriangleIntersection){
            glm::vec3 location = direction * min_intersection;
            if(intersectionTest)
                return std::make_pair(location, 0);
            
            return std::make_pair(barycentric_coordinates, 0);
        } else if(hasSphereIntersection) {
            glm::vec3 location = direction * min_intersection;
            if(intersectionTest)
                return std::make_pair(location, 1);
            
            if(min_intersection < 0.000001) //to avoid reflection intersections due to floating point operation inaccuracy
                return std::make_pair(glm::vec3(-1, -1, -1), 2);
            
            location = initial + direction*min_intersection;
            return std::make_pair(location, 1);
        } else {
           return std::make_pair(glm::vec3(-1, -1, -1), 2);
        }
        
        
    } else
        return std::make_pair(glm::vec3(-1, -1, -1), 2);
    
    
    
}

bool doesIntersect(glm::vec3 direction, glm::vec3 initial, float toLightMagnitude) { //for shadow rays
    //intersection calculations
    Triangle minTriangle;
    Sphere smallestSphere;
    
    
    bool intersectionTest = true;
    std::pair<glm::vec3, int> intersectionDetails =  intersectCalculations(direction, initial, minTriangle, smallestSphere, intersectionTest);


    
    
    if(intersectionDetails.second !=2) {
        float magnitude = glm::length(intersectionDetails.first);
        if(magnitude < (toLightMagnitude))
            return true;
        else
            return false;
    } else
        return false;
    
    
    
}




void sphereIntersecttoPhong(glm::vec3 &phong_position, glm::vec3 &phong_diffuse, glm::vec3 &phong_normal, glm::vec3 &phong_specular, float &phong_shininess, Sphere  & smallestSphere) {
    //calculate normal for sphere and assign to phong normal
    glm::vec3 center_point = glm::vec3(smallestSphere.position[0], smallestSphere.position[1], smallestSphere.position[2]);
    glm::vec3 normalVec = glm::normalize(phong_position - center_point);
    
    
    //phong parameters
    
    phong_normal = normalVec;
    phong_diffuse = glm::vec3(smallestSphere.color_diffuse[0], smallestSphere.color_diffuse[1], smallestSphere.color_diffuse[2]);
    phong_specular = glm::vec3(smallestSphere.color_specular[0], smallestSphere.color_specular[1], smallestSphere.color_specular[2]);
    phong_shininess = smallestSphere.shininess;
    
}

void barycentrictoPhong(glm::vec3 barycentric_coordinates, glm::vec3 &phong_position, glm::vec3 &phong_diffuse, glm::vec3 &phong_normal, glm::vec3 &phong_specular, float &phong_shininess, Triangle  minTriangle){
    glm::vec3 normal_acontrib = barycentric_coordinates.x * glm::vec3(minTriangle.v[0].normal[0], minTriangle.v[0].normal[1], minTriangle.v[0].normal[2]);
    glm::vec3 normal_bcontrib = barycentric_coordinates.y * glm::vec3(minTriangle.v[1].normal[0], minTriangle.v[1].normal[1], minTriangle.v[1].normal[2]);
    glm::vec3 normal_ccontrib = barycentric_coordinates.z * glm::vec3(minTriangle.v[2].normal[0], minTriangle.v[2].normal[1], minTriangle.v[2].normal[2]);
    
    phong_normal = normal_acontrib + normal_bcontrib + normal_ccontrib;
    
    
    glm::vec3 diffuse_acontrib = barycentric_coordinates.x * glm::vec3(minTriangle.v[0].color_diffuse[0], minTriangle.v[0].color_diffuse[1], minTriangle.v[0].color_diffuse[2]);
    glm::vec3 diffuse_bcontrib = barycentric_coordinates.y * glm::vec3(minTriangle.v[1].color_diffuse[0], minTriangle.v[1].color_diffuse[1], minTriangle.v[1].color_diffuse[2]);
    glm::vec3 diffuse_ccontrib = barycentric_coordinates.z * glm::vec3(minTriangle.v[2].color_diffuse[0], minTriangle.v[2].color_diffuse[1], minTriangle.v[2].color_diffuse[2]);
    
    phong_diffuse = diffuse_acontrib + diffuse_bcontrib + diffuse_ccontrib;
    
    glm::vec3 specular_acontrib = barycentric_coordinates.x * glm::vec3(minTriangle.v[0].color_specular[0], minTriangle.v[0].color_specular[1], minTriangle.v[0].color_specular[2]);
    glm::vec3 specular_bcontrib = barycentric_coordinates.y * glm::vec3(minTriangle.v[1].color_specular[0], minTriangle.v[1].color_specular[1], minTriangle.v[1].color_specular[2]);
    glm::vec3 specular_ccontrib = barycentric_coordinates.z * glm::vec3(minTriangle.v[2].color_specular[0], minTriangle.v[2].color_specular[1], minTriangle.v[2].color_specular[2]);
    
    phong_specular = specular_acontrib + specular_bcontrib + specular_ccontrib;
    
    phong_shininess = barycentric_coordinates.x * minTriangle.v[0].shininess + barycentric_coordinates.y * minTriangle.v[1].shininess + barycentric_coordinates.z * minTriangle.v[2].shininess;
    glm::vec3 vertexa = glm::vec3(minTriangle.v[0].position[0], minTriangle.v[0].position[1], minTriangle.v[0].position[2]);
    glm::vec3 vertexb = glm::vec3(minTriangle.v[1].position[0], minTriangle.v[1].position[1], minTriangle.v[1].position[2]);
    glm::vec3 vertexc = glm::vec3(minTriangle.v[2].position[0], minTriangle.v[2].position[1], minTriangle.v[2].position[2]);
    
    phong_position = barycentric_coordinates.x * vertexa + barycentric_coordinates.y * vertexb + barycentric_coordinates.z * vertexc;
}

glm::vec3 raytoColor(glm::vec3 direction, glm::vec3 initial, glm::vec3 prev_specular, bool  & intersectionSuccess) {
    
    float recursedSpecular = glm::length(prev_specular);
    if(recursedSpecular < 0.0001)
        return glm::vec3(0, 0, 0);
    
    //sphere intersection
    Sphere smallestSphere;
    Triangle minTriangle;
    bool intersectionTest = false;
    std::pair<glm::vec3, int> intersectionDetails = intersectCalculations(direction, initial, minTriangle, smallestSphere, intersectionTest);
    
    
    // intersection evaluation calculations
    glm::vec3 phong_normal;
    glm::vec3 phong_diffuse;
    glm::vec3 phong_specular;
    glm::vec3 phong_position;
    float phong_shininess;
    
    
    
    
    //intersection calculations
    
    if(intersectionDetails.second == 0) {
        barycentrictoPhong(intersectionDetails.first, phong_position, phong_diffuse, phong_normal, phong_specular, phong_shininess, minTriangle);
    } else if (intersectionDetails.second == 1) {
        phong_position = intersectionDetails.first;
        sphereIntersecttoPhong(phong_position, phong_diffuse, phong_normal, phong_specular, phong_shininess, smallestSphere);
    }
    
    
    
    
    
    
    
    
    
    glm::vec3 pixelLightColor = glm::vec3(0, 0, 0);
    if(!(intersectionDetails.second == 2)) {
        for(int i=0; i<num_lights; ++i) { //evaluate phong model for pixel color coming from all lights
            
            Light currentLight = lights[i];
            glm::vec3 lightVector = glm::vec3(currentLight.position[0], currentLight.position[1], currentLight.position[2]) - phong_position;
            float magnitude = glm::length(lightVector);
            lightVector = glm::normalize(lightVector);
            
            bool doesIntersectBool = doesIntersect(lightVector, phong_position, magnitude); //shadow ray test
            
            //attenuation of shadows
            float a = 0;
            float b = 0.01;
            float c = 0.05;
            float attenuation = 0;
            if(doesIntersectBool) {
                attenuation = a*1 + b*magnitude + c*magnitude*magnitude;
                attenuation = 1/(attenuation);
            }
            
            
            
            
            
            glm::vec3 currentLightColor = glm::vec3(currentLight.color[0], currentLight.color[1], currentLight.color[2]);
            glm::vec3 reflectionVector = (lightVector * -1.0f) + (glm::dot(lightVector, phong_normal))*(phong_normal * 2.0f);
            
            float diffuseFactor = glm::dot(lightVector, phong_normal);
            if (diffuseFactor < 0)
                diffuseFactor = 0;
            
            glm::vec3 phong_direction = direction * -1.0f;
            
            float specularFactor = glm::dot(reflectionVector, phong_direction);
            if(specularFactor < 0)
                specularFactor = 0;
            
            
            
            float specFactComplete = pow(specularFactor, phong_shininess);
            
            glm::vec3 currentLightValue = currentLightColor * (phong_diffuse * diffuseFactor + phong_specular * specFactComplete);
            currentLightValue -= glm::vec3(attenuation, attenuation, attenuation);
            
            if(currentLightValue.x<0) //less than checks are due to shadows
                currentLightValue.x = 0;
            
            if (currentLightValue.y < 0)
                currentLightValue.y = 0;
            
            if (currentLightValue.z < 0)
                currentLightValue.z = 0;
            
            if(currentLightValue.x>1) //greater than checks due to adding lights
                currentLightValue.x = 1;
            
            if (currentLightValue.y > 1)
                currentLightValue.y = 1;
            
            if (currentLightValue.z > 1)
                currentLightValue.z = 1;
            
            
            
            pixelLightColor += currentLightValue;
            
        }
    }
    
    
    
    
    float r = (pixelLightColor.x) * 255;
    float g = (pixelLightColor.y)  * 255;
    float b = (pixelLightColor.z) * 255;
   
    
    if(!(intersectionDetails.second == 2)) {
        intersectionSuccess = true;
        glm::vec3 reverseDirection = direction * -1.0f;
        glm::vec3 reflectedRay = reverseDirection * -1.0f + (phong_normal)*(glm::dot(reverseDirection, phong_normal))*2.0f;
        glm::vec3 specular_prev = phong_specular * prev_specular;
        bool didIntersect = false;
        
        glm::vec3 colorfromRay = raytoColor(reflectedRay, phong_position, specular_prev, didIntersect);
        
        if(didIntersect) {
            r= (1-phong_specular.x)*(r) + phong_specular.x * colorfromRay.x;
            g= (1-phong_specular.y)*(g) + phong_specular.y * colorfromRay.y;
            b= (1-phong_specular.z)*(b) + phong_specular.z * colorfromRay.z;
        }
    }
    
    return glm::vec3(r, g, b);
}



//MODIFY THIS FUNCTION
void draw_scene()
{
    //shoot ray first parameters
    int z = -1;
    double param, result;
    param = 30.0;
    result = tan (param*(M_PI/180));
    
    
    
    //top right
    float ty = result;
    float aspectRatio = (640/480.0);
    float rx = aspectRatio * result;
    
    
    
    //left
    float lx = -1 * aspectRatio * result;
    float lxOriginal = lx;
    
    //bottom
    float by = -1 * result;
    
    
    float incrementLR = (rx - lx)/640;
    float incrementBT = (ty - by)/480;
    
    
    //a simple test output
    for(unsigned int y=0; y<HEIGHT; y++)
    {
        
        
        
        glPointSize(2.0);
        glBegin(GL_POINTS);
        
        
        for(unsigned int x=0; x<WIDTH; x++)
        {
            
            
            
            
            
            
            int r = 0;
            int g = 0;
            int b = 0;
            
            glm::vec3 color = glm::vec3(0, 0, 0);
            for(int i=0; i<5; ++i) { //pixel split into four for and averaged (supersampling for antialiasing)
                for(int j=0; j<5; ++j) {
                    float xaddition = -0.5 + i*0.25;
                    float yaddition = -0.5 + j*0.25;
                    
                    
                    //for ray
                    glm::vec3 ray = glm::vec3(lx + xaddition*incrementLR, by + yaddition*incrementBT, z);
                    glm::vec3 direction = glm::normalize(ray);
                    
                    bool didIntersect = false;
                    color += raytoColor(direction , glm::vec3(0,0,0), glm::vec3(1, 1, 1), didIntersect);
                    
                    
                }
            }
            
            lx += incrementLR;
            color = color/25.0f;
            //add ambient light to each pixel
            color.x +=  (ambient_light[0]*255);
            color.y += (ambient_light[1] *255);
            color.z += (ambient_light[2] * 255);
            
            
            if(color.x>255) //greater than checks due to adding ambient light
                color.x = 255;
            
            if(color.y>255)
                color.y = 255;
            
            if(color.z>255)
                color.z = 255;
            
            r = roundf(color.x);
            g = roundf(color.y);
            b = roundf(color.z);
            
            unsigned char rchar = r;
            unsigned char gchar = g;
            unsigned char bchar = b;
            plot_pixel(x, y, rchar, gchar, bchar);
            
        }
        
        glEnd();
        glFlush();
        by +=incrementBT;
        lx = lxOriginal;
        
    }
    
    
    
    
    printf("Done!\n");
    fflush(stdout);
}





void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
    glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
    glVertex2i(x,y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
    buffer[y][x][0] = r;
    buffer[y][x][1] = g;
    buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
    plot_pixel_display(x,y,r,g,b);
    if(mode == MODE_JPEG)
        plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
    printf("Saving JPEG file: %s\n", filename);
    
    ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
    if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
        printf("Error in Saving\n");
    else
        printf("File saved Successfully\n");
}

void parse_check(const char *expected, char *found)
{
    if(strcasecmp(expected,found))
    {
        printf("Expected '%s ' found '%s '\n", expected, found);
        printf("Parse error, abnormal abortion\n");
        exit(0);
    }
}

void parse_doubles(FILE* file, const char *check, double p[3])
{
    char str[100];
    fscanf(file,"%s",str);
    parse_check(check,str);
    fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
    printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE *file, double *r)
{
    char str[100];
    fscanf(file,"%s",str);
    parse_check("rad:",str);
    fscanf(file,"%lf",r);
    printf("rad: %f\n",*r);
}

void parse_shi(FILE *file, double *shi)
{
    char s[100];
    fscanf(file,"%s",s);
    parse_check("shi:",s);
    fscanf(file,"%lf",shi);
    printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
    FILE * file = fopen(argv,"r");
    int number_of_objects;
    char type[50];
    Triangle t;
    Sphere s;
    Light l;
    fscanf(file,"%i", &number_of_objects);
    
    printf("number of objects: %i\n",number_of_objects);
    
    parse_doubles(file,"amb:",ambient_light);
    
    for(int i=0; i<number_of_objects; i++)
    {
        fscanf(file,"%s\n",type);
        printf("%s\n",type);
        if(strcasecmp(type,"triangle")==0)
        {
            printf("found triangle\n");
            for(int j=0;j < 3;j++)
            {
                parse_doubles(file,"pos:",t.v[j].position);
                parse_doubles(file,"nor:",t.v[j].normal);
                parse_doubles(file,"dif:",t.v[j].color_diffuse);
                parse_doubles(file,"spe:",t.v[j].color_specular);
                parse_shi(file,&t.v[j].shininess);
            }
            
            if(num_triangles == MAX_TRIANGLES)
            {
                printf("too many triangles, you should increase MAX_TRIANGLES!\n");
                exit(0);
            }
            triangles[num_triangles++] = t;
        }
        else if(strcasecmp(type,"sphere")==0)
        {
            printf("found sphere\n");
            
            parse_doubles(file,"pos:",s.position);
            parse_rad(file,&s.radius);
            parse_doubles(file,"dif:",s.color_diffuse);
            parse_doubles(file,"spe:",s.color_specular);
            parse_shi(file,&s.shininess);
            
            if(num_spheres == MAX_SPHERES)
            {
                printf("too many spheres, you should increase MAX_SPHERES!\n");
                exit(0);
            }
            spheres[num_spheres++] = s;
        }
        else if(strcasecmp(type,"light")==0)
        {
            printf("found light\n");
            parse_doubles(file,"pos:",l.position);
            parse_doubles(file,"col:",l.color);
            
            if(num_lights == MAX_LIGHTS)
            {
                printf("too many lights, you should increase MAX_LIGHTS!\n");
                exit(0);
            }
            lights[num_lights++] = l;
        }
        else
        {
            printf("unknown type in scene description:\n%s\n",type);
            exit(0);
        }
    }
    return 0;
}

void display()
{
}

void init()
{
    glMatrixMode(GL_PROJECTION);
    glOrtho(0,WIDTH,0,HEIGHT,1,-1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    glClearColor(0,0,0,0);
    glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
    //hack to make it only draw once
    static int once=0;
    if(!once)
    {
        draw_scene();
        if(mode == MODE_JPEG)
            save_jpg();
    }
    once=1;
}

int main(int argc, char ** argv)
{
    if ((argc < 2) || (argc > 3))
    {
        printf ("Usage: %s <input scenefile> [output jpegname]\n", argv[0]);
        exit(0);
    }
    if(argc == 3)
    {
        mode = MODE_JPEG;
        filename = argv[2];
    }
    else if(argc == 2)
        mode = MODE_DISPLAY;
    
    glutInit(&argc,argv);
    loadScene(argv[1]);
    
    glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
    glutInitWindowPosition(0,0);
    glutInitWindowSize(WIDTH,HEIGHT);
    int window = glutCreateWindow("Ray Tracer");
    glutDisplayFunc(display);
    glutIdleFunc(idle);
    init();
    glutMainLoop();
}

