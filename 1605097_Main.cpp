#include<windows.h>
#include <GL/glut.h>
#include <bits/stdc++.h>
#include "1605097_header.h"
#include "bitmap_image.hpp"
#define pi (2*acos(0.0))
#define window_width 500
#define window_height 500
#define view_angle 80
using namespace std;
double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle;
int x_as=0;
int y_as=0;
int z_as=0;
point position,u,r,l;
vector<string> separated_strings[500];
int width_image,height_image;
int recursion_level;
void drawAxes()
{

    glColor3f(1.0, 1.0, 1.0);
    glBegin(GL_LINES);
    {
        glVertex3f( 100,0,0);
        glVertex3f(-100,0,0);

        glVertex3f(0,-100,0);
        glVertex3f(0, 100,0);

        glVertex3f(0,0, 100);
        glVertex3f(0,0,-100);
    }
    glEnd();

}


void drawGrid()
{
    int i;
    if(drawgrid==1)
    {
        glColor3f(0.6, 0.6, 0.6);	//grey
        glBegin(GL_LINES);
        {
            for(i=-8; i<=8; i++)
            {

                if(i==0)
                    continue;	//SKIP the MAIN axes

                //lines parallel to Y-axis
                glVertex3f(i*10, -90, 0);
                glVertex3f(i*10,  90, 0);

                //lines parallel to X-axis
                glVertex3f(-90, i*10, 0);
                glVertex3f( 90, i*10, 0);
            }
        }
        glEnd();
    }
}





void drawLine(double a)
{
    glColor3f(0,128,0);
    glBegin(GL_LINES);
    {
        glVertex3f(-a,a,0);
        glVertex3f(a,a,0);
        glVertex3f(-a,-a,0);
        glVertex3f(-a,a,0);
        glVertex3f(-a,-a,0);
        glVertex3f(a,-a,0);
        glVertex3f(a,-a,0);
        glVertex3f(a,a,0);
    }
    glEnd();
}
void Normalize(struct point vector1)
{

    double mag=sqrt(vector1.x*vector1.x+vector1.y*vector1.y+vector1.z*vector1.z);
    vector1.x=vector1.x/mag;
    vector1.y=vector1.y/mag;
    vector1.z=vector1.z/mag;

}
void rotate_left()
{
    double rangle;
    rangle=(angle*1.0*pi)/180;
    point cross_mult_1,cross_mult_2;
    cross_mult_1=CrossProduct(u,l);
    //printf(" vector is %f %f %f",cross_mult_1.x,cross_mult_1.y,cross_mult_1.z);
    cross_mult_2=CrossProduct(u,r);
    //printf(" vector is %f %f %f",l.x,l.y,l.z);
    l.x=l.x*cos(rangle)+cross_mult_1.x*sin(rangle);
    l.y=l.y*cos(rangle)+cross_mult_1.y*sin(rangle);
    l.z=l.z*cos(rangle)+cross_mult_1.z*sin(rangle);
    r.x=r.x*cos(rangle)+cross_mult_2.x*sin(rangle);
    r.y=r.y*cos(rangle)+cross_mult_2.y*sin(rangle);
    r.z=r.z*cos(rangle)+cross_mult_2.z*sin(rangle);
    Normalize(l);
    //printf("After normalizing vector is %f %f %f",l.x,l.y,l.z);
    Normalize(r);



}
void rotate_right()
{
    double rangle;
    rangle=(angle*1.0*pi)/180;
    point cross_mult_1,cross_mult_2;
    cross_mult_1=CrossProduct(u,l);
    cross_mult_2=CrossProduct(u,r);
    l.x=l.x*cos(rangle)-cross_mult_1.x*sin(rangle);
    l.y=l.y*cos(rangle)-cross_mult_1.y*sin(rangle);
    l.z=l.z*cos(rangle)-cross_mult_1.z*sin(rangle);
    r.x=r.x*cos(rangle)-cross_mult_2.x*sin(rangle);
    r.y=r.y*cos(rangle)-cross_mult_2.y*sin(rangle);
    r.z=r.z*cos(rangle)-cross_mult_2.z*sin(rangle);
    Normalize(l);
    Normalize(r);



}
void look_down()
{
    double rangle;
    rangle=(angle*1.0*pi)/180;
    point cross_mult_1,cross_mult_2;
    cross_mult_1=CrossProduct(r,l);
    cross_mult_2=CrossProduct(r,u);
    l.x=l.x*cos(rangle)-cross_mult_1.x*sin(rangle);
    l.y=l.y*cos(rangle)-cross_mult_1.y*sin(rangle);
    l.z=l.z*cos(rangle)-cross_mult_1.z*sin(rangle);
    u.x=u.x*cos(rangle)-cross_mult_2.x*sin(rangle);
    u.y=u.y*cos(rangle)-cross_mult_2.y*sin(rangle);
    u.z=u.z*cos(rangle)-cross_mult_2.z*sin(rangle);
    Normalize(l);
    Normalize(u);

}
void look_up()
{
    double rangle;
    rangle=(angle*1.0*pi)/180;
    point cross_mult_1,cross_mult_2;
    cross_mult_1=CrossProduct(r,l);
    cross_mult_2=CrossProduct(r,u);
    l.x=l.x*cos(rangle)+cross_mult_1.x*sin(rangle);
    l.y=l.y*cos(rangle)+cross_mult_1.y*sin(rangle);
    l.z=l.z*cos(rangle)+cross_mult_1.z*sin(rangle);
    u.x=u.x*cos(rangle)+cross_mult_2.x*sin(rangle);
    u.y=u.y*cos(rangle)+cross_mult_2.y*sin(rangle);
    u.z=u.z*cos(rangle)+cross_mult_2.z*sin(rangle);
    Normalize(l);
    Normalize(u);

}
void anti_clock()
{
    double rangle;
    rangle=(angle*1.0*pi)/180;
    point cross_mult_1,cross_mult_2;
    cross_mult_1=CrossProduct(l,u);
    cross_mult_2=CrossProduct(l,r);
    u.x=u.x*cos(rangle)+cross_mult_1.x*sin(rangle);
    u.y=u.y*cos(rangle)+cross_mult_1.y*sin(rangle);
    u.z=u.z*cos(rangle)+cross_mult_1.z*sin(rangle);
    r.x=r.x*cos(rangle)+cross_mult_2.x*sin(rangle);
    r.y=r.y*cos(rangle)+cross_mult_2.y*sin(rangle);
    r.z=r.z*cos(rangle)+cross_mult_2.z*sin(rangle);
    Normalize(l);
    Normalize(u);
}
void clockR()
{
    double rangle;
    rangle=(angle*1.0*pi)/180;
    point cross_mult_1,cross_mult_2;
    cross_mult_1=CrossProduct(l,u);
    cross_mult_2=CrossProduct(l,r);
    u.x=u.x*cos(rangle)-cross_mult_1.x*sin(rangle);
    u.y=u.y*cos(rangle)-cross_mult_1.y*sin(rangle);
    u.z=u.z*cos(rangle)-cross_mult_1.z*sin(rangle);
    r.x=r.x*cos(rangle)-cross_mult_2.x*sin(rangle);
    r.y=r.y*cos(rangle)-cross_mult_2.y*sin(rangle);
    r.z=r.z*cos(rangle)-cross_mult_2.z*sin(rangle);
    Normalize(l);
    Normalize(u);
}

void specialKeyListener(int key, int x,int y)
{
    float amount= 3.0f;
    switch(key)
    {
    case GLUT_KEY_DOWN:		//down arrow key
        position.x=position.x-l.x*amount;
        position.y=position.y-l.y*amount;
        break;
    case GLUT_KEY_UP:		// up arrow key
        position.x=position.x+l.x*amount;
        position.y=position.y+l.y*amount;
        break;

    case GLUT_KEY_RIGHT:
        position.x=position.x+r.x*amount;
        position.y=position.y+r.y*amount;
        break;
    case GLUT_KEY_LEFT:
        position.x=position.x-r.x*amount;
        position.y=position.y-r.y*amount;
        break;

    case GLUT_KEY_PAGE_UP:
        position.z=position.z+u.z*3.0;
        // check_z=position.z;
        break;
    case GLUT_KEY_PAGE_DOWN:
        position.z=position.z-u.z*3.0;
        // check_z=position.z;
        break;

    case GLUT_KEY_INSERT:
        break;

    case GLUT_KEY_HOME:
        break;
    case GLUT_KEY_END:
        break;

    default:
        break;
    }
}


void mouseListener(int button, int state, int x, int y) 	//x, y is the x-y of the screen (2D)
{
    switch(button)
    {
    case GLUT_LEFT_BUTTON:
        if(state == GLUT_DOWN) 		// 2 times?? in ONE click? -- solution is checking DOWN or UP
        {
            drawaxes=1-drawaxes;
        }
        break;

    case GLUT_RIGHT_BUTTON:
        //........
        break;

    case GLUT_MIDDLE_BUTTON:
        //........
        break;

    default:
        break;
    }
}
void remove_mem()
{
    for(int i=0 ;i<object.size(); i++)
    {
        delete object[i];
    }
    object.clear();
    for(int i=0;i<light.size();i++){
        delete light[i];
    }
    light.clear();

}
void loadData()
{
    //cout<<"up to date.........."<<endl;
    string type;
    ifstream scene;
    scene.open("G:\\Academic_Things\\4-1\\Graphics Sessional-410\\Offline-3\\Ray_Tracing\\scene.txt");
    if (!scene.is_open())
    {
        cout << "Error while opening" << endl;
        exit(1);
    }
    string lines;
    int count_line = 0;
    while (getline(scene, lines))
    {
        stringstream s(lines);
        string tokens;
        while (s >> tokens)
        {
            //cout<<tokens<<"...tokens..."<<count_line<<"...countline..."<<endl;
            separated_strings[count_line].push_back(tokens);

        }

        count_line++;
    }
    double d=stod(separated_strings[50][0].c_str());
    cout<<separated_strings[50][0]<<".........hi"<<d<<endl;
    Object *temp;
    Light *temp_light;
    temp=new checkerBoard(1000,20);
    temp->setCoefficients(0.4,0.4,0.4,0.1);
    temp->setShine(1.0);
    object.push_back(temp);
    int total_objects;

    for (int i = 0; i < 5; i++)
    {
        if (i == 0)
        {
            recursion_level=stod(separated_strings[i][0].c_str());
        }
        if(i==1)
        {
            width_image=stod(separated_strings[i][0].c_str());
            height_image=stod(separated_strings[i][0].c_str());
        }
        if(i==3)
        {
            total_objects=stod(separated_strings[i][0].c_str());
        }

    }
    string command_line;
    int light_no;
    for (int i=4; i<count_line; i++)
    {
        for (int j = 0; j < separated_strings[i].size(); j++)
        {
            command_line = separated_strings[i][j];
            if(command_line=="general")
            {
                //cout<<command_line<<"..."<<i+1<<endl;
                double a,b,c,d,e,f,g,h,io,j;
                point ref_point;
                double length,width,height;
                a=stod(separated_strings[i + 1][0].c_str());
                b=stod(separated_strings[i + 1][1].c_str());
                c=stod(separated_strings[i + 1][2].c_str());
                d=stod(separated_strings[i + 1][3].c_str());
                e=stod(separated_strings[i + 1][4].c_str());
                f=stod(separated_strings[i + 1][5].c_str());
                g=stod(separated_strings[i + 1][6].c_str());
                h=stod(separated_strings[i + 1][7].c_str());
                io=stod(separated_strings[i + 1][8].c_str());
                j=stod(separated_strings[i + 1][9].c_str());
                ref_point.x=stod(separated_strings[i + 2][0].c_str());
                ref_point.y=stod(separated_strings[i + 2][1].c_str());
                ref_point.z=stod(separated_strings[i + 2][2].c_str());
                length=stod(separated_strings[i + 2][3].c_str());
                width=stod(separated_strings[i + 2][4].c_str());
                height=stod(separated_strings[i + 2][5].c_str());
                double red=stod(separated_strings[i + 3][0].c_str());
                double green=stod(separated_strings[i + 3][1].c_str());
                double blue=stod(separated_strings[i + 3][2].c_str());
                temp=new General(a,b,c,d,e,f,g,h,io,j,ref_point,length,width,height);
                temp->setColour(red,green,blue);
                double ambient, diffuse, reflection, specular, shininess;
                ambient=stod(separated_strings[i + 4][0].c_str());
                diffuse=stod(separated_strings[i + 4][1].c_str());
                reflection=stod(separated_strings[i + 4][2].c_str());
                specular=stod(separated_strings[i + 4][3].c_str());
                shininess=stod(separated_strings[i + 5][0].c_str());
                temp->setCoefficients(ambient,diffuse,reflection,specular);
                temp->setShine(shininess);
                object.push_back(temp);


            }
            else if(command_line=="sphere")
            {
                cout<<command_line<<"..."<<i<<endl;
                double center_X,center_Y,center_Z,radius;
                center_X=stod(separated_strings[i + 1][0].c_str());
                center_Y=stod(separated_strings[i + 1][1].c_str());
                center_Z=stod(separated_strings[i + 1][2].c_str());
                radius=stod(separated_strings[i + 2][0].c_str());
                point ref_point(center_X,center_Y,center_Z);
                temp=new Sphere(ref_point,radius);
                double r,g,b;
                r=stod(separated_strings[i + 3][0].c_str());
                g=stod(separated_strings[i + 3][1].c_str());
                b=stod(separated_strings[i + 3][2].c_str());
                temp->setColour(r,g,b);
                double ambient, diffuse, reflection, specular, shininess;
                ambient=stod(separated_strings[i + 4][0].c_str());
                diffuse=stod(separated_strings[i + 4][1].c_str());
                reflection=stod(separated_strings[i + 4][2].c_str());
                specular=stod(separated_strings[i + 4][3].c_str());
                shininess=stod(separated_strings[i + 5][0].c_str());
                temp->setCoefficients(ambient,diffuse,reflection,specular);
                temp->setShine(shininess);
              //  temp->print();
                object.push_back(temp);

            }
            else if(command_line=="triangle")
            {
                cout<<command_line<<"..."<<i<<endl;
                point a,b,c;
                a.x=stod(separated_strings[i + 1][0].c_str());
                a.y=stod(separated_strings[i + 1][1].c_str());
                a.z=stod(separated_strings[i + 1][2].c_str());
                b.x=stod(separated_strings[i + 2][0].c_str());
                b.y=stod(separated_strings[i + 2][1].c_str());
                b.z=stod(separated_strings[i + 2][2].c_str());
                c.x=stod(separated_strings[i + 3][0].c_str());
                c.y=stod(separated_strings[i + 3][1].c_str());
                c.z=stod(separated_strings[i + 3][2].c_str());
                temp=new Triangle(a,b,c);
                double r,g,blue;
                r=stod(separated_strings[i + 4][0].c_str());
                g=stod(separated_strings[i + 4][1].c_str());
                blue=stod(separated_strings[i + 4][2].c_str());
                temp->setColour(r,g,blue);
                double ambient, diffuse, reflection, specular, shininess;
                ambient=stod(separated_strings[i + 5][0].c_str());
                diffuse=stod(separated_strings[i + 5][1].c_str());
                reflection=stod(separated_strings[i + 5][2].c_str());
                specular=stod(separated_strings[i + 5][3].c_str());
                shininess=stod(separated_strings[i + 6][0].c_str());
                temp->setCoefficients(ambient,diffuse,reflection,specular);
                temp->setShine(shininess);
                temp->print();
                object.push_back(temp);

            }



        }



    }
    cout<<"light_no:..."<<endl;
    cin>>light_no;
    for(int i=0; i<2*light_no; i=i+2)
    {
        point l;
        l.x=stod(separated_strings[count_line-2*light_no+i][0].c_str());
        l.y=stod(separated_strings[count_line-2*light_no+i][1].c_str());
        l.z=stod(separated_strings[count_line-2*light_no+i][2].c_str());
        double r,g,b;
        r=stod(separated_strings[count_line-2*light_no+i+1][0].c_str());
        g=stod(separated_strings[count_line-2*light_no+i+1][1].c_str());
        b=stod(separated_strings[count_line-2*light_no+i+1][2].c_str());
        //cout<<r<<"...r.."<<g<<"...g..."<<b<<"....b..."<<endl;
        temp_light=new Light(l,r,g,b);
      //  temp_light->print();
        light.push_back(temp_light);



    }


}
void capture()
{
    point** img_buffer=new point*[width_image];
    for(int i=0; i<width_image; i++)
    {
        img_buffer[i]=new point[width_image];
    }
    for(int i=0; i<width_image; i++)
    {
        for(int j=0; j<width_image; j++)
        {
            img_buffer[i][j]=point(0,0,0);
        }
    }
    double planeDistance=(window_height/2.0)/(tan((view_angle*pi)/360));
    point temp1=l*planeDistance;
    point top_left=position+l*planeDistance-r*(window_width/2)+u*(window_height/2);
    double du=(window_width*1.0/width_image);
    double dv=(window_height*1.0/height_image);
    top_left=top_left+r*(0.5*du)-u*(0.5*dv);
   // top_left.printPoint();
    int nearest;
    bitmap_image image(width_image,height_image);
    for(int i=0; i<width_image; i++)
    {
        for(int j=0l; j<height_image; j++)
        {
            image.set_pixel(i,j,0,0,0);
        }
    }
    for(int i=0; i<width_image; i++)
    {
        for(int j=0; j<height_image; j++)
        {
            int ri=i*du;
            int ui=j*dv;
            point dir_top_left=top_left+r*ri-u*ui;
            point direction=dir_top_left-position;
            Ray ray(position,direction);
            double *color=new double[3];
            int nearest=-1;
            double min_t=10000000;
            for(int k=0; k<object.size(); k++)
            {
                double t = object[k]->intersection(&ray,color,0);

                if(t>0 &&t<min_t )
                {
                  //  cout<<color[1]<<"..."<<color[2]<<"..."<<color[3]<<endl;

                    min_t = t;
                    nearest=k;
                }

            }
            if(nearest!=-1)
            {

              //cout<<"nearest...."<<endl;

                min_t = object[nearest]->intersection(&ray,color,3);
                if(color[0]<0.0){
                    color[0]=0.0;
                }
                if(color[0]>1.0){
                    color[0]=1.0;
                }
                 if(color[1]<0.0){
                    color[1]=0.0;
                }
                if(color[1]>1.0){
                    color[1]=1.0;
                }
                 if(color[2]<0.0){
                    color[2]=0.0;
                }
                if(color[2]>1.0){
                    color[2]=1.0;
                }
                image.set_pixel(i,j,color[0]*255, color[1]*255,color[2]*255);

            }
        //    delete color;


        }
    }

    cout<<"dhukseee...."<<endl;
    image.save_image("G:\\Academic_Things\\4-1\\Graphics Sessional-410\\Offline-3\\Ray_Tracing\\output.bmp");



}
void keyboardListener(unsigned char key, int x,int y)
{
    switch(key)
    {
    case '0':
        capture();

    case '1':
        rotate_left();
        break;
    case '2':
        rotate_right();
        break;
    case '3':
        look_up();
        break;
    case '4':
        look_down();
        break;
    case '5':
        clockR();
        break;
    case '6':
        anti_clock();
        break;
    default:
        break;
    }
}

void display()
{

    //clear the display
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0,0,0,0);	//color black
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(	position.x, position.y, position.z,
                position.x+l.x, position.y+l.y,  position.z+l.z,
                u.x, u.y, u.z);

    glMatrixMode(GL_MODELVIEW);
    drawAxes();
    drawGrid();
    for(int i = 0; i < object.size(); i++)
    {

        object[i]->draw();
    }
    for(int i = 0; i < light.size(); i++)
    {

        light[i]->draw();
    }
    glutSwapBuffers();
}


void animate()
{
    glutPostRedisplay();
}

void init()
{

    loadData();
    atexit(remove_mem);
    drawgrid=0;
    drawaxes=1;
    position.x=100;
    position.y=100;
    position.z=0;
    u.x=0;
    u.y=0;
    u.z=1;
    r.x=-1/sqrt(2);
    r.y=1/sqrt(2);
    r.z=0;
    l.x=-1/sqrt(2);
    l.y=-1/sqrt(2);
    l.z=0;
    printf("init  e vector is %f %f %f",l.x,l.y,l.z);
    angle = 3.0f;
    glClearColor(0,0,0,0);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(80,	1,	1,	1000.0);



}

int main(int argc, char **argv)
{
    glutInit(&argc,argv);
    glutInitWindowSize(500, 500);
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

    glutCreateWindow("My OpenGL Program");

    init();
    //capture();
    glEnable(GL_DEPTH_TEST);	//enable Depth Testing

    glutDisplayFunc(display);	//display callback function
    glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring
    //glutSpecialFunc(specialKeyListener);
    glutKeyboardFunc(keyboardListener);
    // printf("Main e vector is %d %d %d",l.x,l.y,l.z);
    glutSpecialFunc(specialKeyListener);
    // printf("Main e vector is %d %d %d",&l.x,&l.y,&l.z);
    glutMouseFunc(mouseListener);

    glutMainLoop();
    //The main loop of OpenGL

    return 0;
}
