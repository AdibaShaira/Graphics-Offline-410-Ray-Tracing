#ifndef BASE_HPP_INCLUDED
#define BASE_HPP_INCLUDED
#define pi (2*acos(0.0))
extern int recursion_level;
int jakhushi=0;
using namespace std;
struct point
{
    double x,y,z;
    point()
    {
        x=0;
        y=0;
        x=0;
    }
    point(double a,double b,double c)
    {
        this->x=a;
        this->y=b;
        this->z=c;
    }
    point operator-(point p)
    {
        point temp;
        temp.x=this->x - p.x;
        temp.y=this->y - p.y;
        temp.z=this->z - p.z;
        return temp;
    }
    point operator+(point p)
    {
        point temp;
        temp.x=this->x + p.x;
        temp.y=this->y + p.y;
        temp.z=this->z + p.z;
        return temp;
    }
    point operator*(double sc)
    {
        point temp;
        temp.x=this->x*sc;
        temp.y=this->y*sc;
        temp.z=this->z*sc;
        return temp;

    }
    point operator/(double sc)
    {
        point temp;
        temp.x=this->x/sc;
        temp.y=this->y/sc;
        temp.z=this->z/sc;
        return temp;

    }
    void printPoint()
    {
        cout<<"x-> "<<x<<"y-> "<<y<<"z-> "<<z<<endl;
    }
    point Normalize(point v)
    {
        double r;
        r = sqrt((v.x * v.x) + (v.y * v.y) + (v.z * v.z));
        v.x = v.x / r;
        v.y = v.y / r;
        v.z = v.z / r;
        return v;
    }
    double Dot_Multiplication(point v2)
    {
        double result;
        result = this->x * v2.x + this->y * v2.y + this->z * v2.z;
        return result;
    }
};
point CrossProduct(struct point vector1, struct point vector2)
{
    struct point result;
    result.x=vector1.y*vector2.z-vector1.z*vector2.y;
    result.y=vector1.z*vector2.x-vector1.x*vector2.z;
    result.z=vector1.x*vector2.y-vector1.y*vector2.x;
    return result;
}
struct color
{
    double r,g,b;
    color()
    {
        r=0;
        g=0;
        b=0;
    }
    color(double R,double G,double B)
    {
        this->r=R;
        this->g=G;
        this->b=B;
    }
    void printColor()
    {
        cout << "r->" << r << " g-> " << g << " b-> " << b<<endl;
    }

};

class Ray
{
public:
    point initial;
    point direction;
    Ray(point i,point d)
    {
        this->initial=i;
        this->direction=d.Normalize(d);
    }
};
class Object
{
public:
    point reference_point;
    double height,width,length;
    double color[3];
    double coEfficient[4];
    int shine;
    Object()
    {
    }
    virtual void draw()
    {
    }
    void setColour(double r,double g,double b)
    {
        color[0]=r;
        color[1]=g;
        color[2]=b;
    }
    void setShine(int shine)
    {
        this->shine=shine;
    }
    void setCoefficients(double a,double b,double c,double d)
    {
        coEfficient[0]=a;
        coEfficient[1]=b;
        coEfficient[2]=c;
        coEfficient[3]=d;
    }
    virtual void print()
    {
    }
    virtual double getD(point p) = 0;
    virtual point getNormal(point p,Ray *r)=0;
    virtual double getintersection(Ray *ray)
    {
        return 0;
    }
    virtual double intersection(Ray *r, double *color, int level)
    {
        return 0;
    }
    point getReflection(Ray *r,point normal)
    {
        double temp=2.0*r->direction.Dot_Multiplication(normal);
        point tempp=normal*temp;
        point reflection=r->direction-tempp;
        point norm=reflection.Normalize(reflection);
        return norm;
    }


};
class Light
{
public:
    point ref_point;
    double colors[3];
    Light(point refp,double r,double g,double b)
    {
        this->ref_point=refp;
        colors[0]=r;
        colors[1]=g;
        colors[2]=b;
    }
    void draw()
    {
        glColor3f(colors[0],colors[1],colors[2]);

        glBegin(GL_QUADS);
        {
            glVertex3f(ref_point.x+4,ref_point.y,ref_point.z+4);
            glVertex3f(ref_point.x+4,ref_point.y,ref_point.z-4);
            glVertex3f(ref_point.x-4,ref_point.y,ref_point.z-4);
            glVertex3f(ref_point.x-4,ref_point.y,ref_point.z+4);
        }
        glEnd();
    }
    void print()
    {
        ref_point.printPoint();
    }
    double getNormal(point p)
    {

    }

};

vector<Object *> object;
vector<Light *> light;
class Sphere:public Object
{
public:
    point center;
    double radius;
    struct point points[100][100];
    int slices, stacks;
    Sphere(point center,double radius)
    {
        this->center=center;
        this->radius=radius;

    }

    void print()
    {

        center.printPoint();
        cout<<"radius...."<<radius<<endl;
        cout<<"colors.."<<color[0]<<"..."<<color[1]<<"..."<<color[2]<<endl;
    }
    void draw()
    {
        glColor3f(color[0], color[1], color[2]);
        glPushMatrix();
        glTranslatef(center.x, center.y, center.z);
        glutSolidSphere(radius, 100, 100);
        glPopMatrix();

    }
    double getD(point p)
    {

    }
    point getNormal(point intersection,Ray *r)
    {
        point normal=intersection-this->center;
        normal=normal.Normalize(normal);
        double ao=(r->initial.x-this->center.x);
        double bo=(r->initial.y-this->center.y);
        double co=(r->initial.z-this->center.z);
        double dis=sqrt(ao*ao+bo*bo+co*co);
        if(dis<this->radius)
        {
            normal=normal*-1.0;
        }
        return normal;
    }
    double getintersection(Ray *ray)
    {
        point ac=ray->initial-this->center;
        double a=1.0;
        double b=2.0*ray->direction.Dot_Multiplication(ac);
        double c=ac.Dot_Multiplication(ac)-(radius*radius);
        //    cout<<"...a..."<<a<<"...b..."<<b<<"...c..."<<c<<endl;

        double dis=pow(b,2)-(4*a*c);

        if(dis<0.0)
        {
            return -1;

        }
        else
        {
            double numerator2 = (-b - sqrt(dis))/(2.0*a);
            double numerator1 = (-b + sqrt(dis))/(2.0*a);
            if(numerator1>=0 && numerator2>=0)
            {
                return min(numerator1,numerator2);
            }
            else
            {
                if(numerator1>=0 && numerator2<0)
                {
                    return numerator1;
                }

                else if(numerator1<0 && numerator2>=0)
                {
                    return numerator2;
                }
                else if(numerator1<0 && numerator2<0)
                {
                    return -1;
                }
            }

        }
    }
    double intersection(Ray *r, double *new_color, int level)
    {

        double t=getintersection(r);
        // cout<<"sphere er t..."<<t<<endl;
        if(t<=0)
        {
            return -1;
        }
        if(level==0)
        {
            return t;
        }
        for(int i=0; i<3; i++)
        {
            new_color[i]=this->color[i]*this->coEfficient[0];
            if(new_color[i]<0.0)
            {
                new_color[i]=0.0;
            }
            else if(new_color[i]>1.0)
            {
                new_color[i]=1.0;
            }

        }

        point temp=r->direction*t;
        point intersectionPoint=r->initial+temp;
        point normal=getNormal(intersectionPoint,r);
        point reflection=getReflection(r,normal);
        double min_t=1000000;
        for(int i=0; i<light.size(); i++)
        {
            point dir=light[i]->ref_point-intersectionPoint;
            point norm_dir=dir.Normalize(dir);
            point initial=intersectionPoint+norm_dir*0.000001;
            Ray ro(initial,norm_dir);
            bool flag=false;
            double t_min = 10000000;
            for(int j=0; j<object.size(); j++)
            {
                double t_temp=object[j]->getintersection(&ro);
                if(t_temp>0 && t_min>t_temp)
                {
                    t_min=t_temp;
                    flag=true;
                    break;

                }
            }
            if(!flag)
            {
                //lambert er jonne
                point d=ro.direction*-1.0;
                double dotn=d.Dot_Multiplication(normal)*2.0;
                point R=d-normal*dotn;
                R=R.Normalize(R);
                point V= r->direction *-1.0;
                V=V.Normalize(V);
                double lambert=max(0.0,ro.direction.Dot_Multiplication(normal));
                double phong = max(0.0,R.Dot_Multiplication(V));
                phong=pow(phong,shine);
                for (int k=0; k<3; k++)
                {
                    new_color[k] =new_color[k]+ light[i]->colors[k]* lambert * coEfficient[1] *this->color[k];
                    new_color[k]=new_color[k]+light[i]->colors[k] * phong * coEfficient[2] *this->color[k];
                    if(new_color[k]<0.0)
                    {
                        new_color[k]=0.0;
                    }
                    else if(new_color[k]>1.0)
                    {
                        new_color[k]=1.0;
                    }

                }

            }

        }
        if(level<recursion_level)
        {
            point initial=intersectionPoint+reflection*0.000001;
            Ray ref_ray(initial,reflection);
            int nearest=-1;
            double *ref_color=new double[3];

            for(int k=0; k<object.size(); k++)
            {
                double to = object[k]->getintersection(&ref_ray);
                if(to>0 && to<min_t )
                {
                    min_t =to;
                    nearest=k;
                }
            }
            if(nearest!=-1)
            {
                min_t = object[nearest]->intersection(&ref_ray,ref_color,level+1);
                for (int k=0; k<3; k++)
                {
                    new_color[k]+=ref_color[k]*coEfficient[3];
                     if(new_color[k]<0.0)
                    {
                        new_color[k]=0.0;
                    }
                    else if(new_color[k]>1.0)
                    {
                        new_color[k]=1.0;
                    }
                }
            }

        }
        return t;



    }


};
class Triangle:public Object
{
public:
    double a,b,c,d;
    point point1;
    point point2;
    point point3;
    Triangle(point p1,point p2,point p3)
    {
        this->point1=p1;
        this->point2=p2;
        this->point3=p3;
        point temp1= this->point2- this->point1;
        point temp2= this->point3- this->point1;
        point cross_mul= CrossProduct(temp1,temp2);
        this->a =  cross_mul.x;
        this->b =  cross_mul.y;
        this->c =  cross_mul.z;
        this->d=cross_mul.x*this->point1.x+cross_mul.y*this->point1.y+cross_mul.z*this->point1.z;

    }
    double getD(point p)
    {
        point temp1,temp2,cross_mul;
        temp1= this->point2- this->point1;
        temp2= this->point3- this->point1;
        cross_mul= CrossProduct(temp1,temp2);
        this->a = cross_mul.x;
        this->b = cross_mul.y;
        this->c = cross_mul.z;
        this->d=cross_mul.x*this->point1.x+cross_mul.y*this->point2.y+cross_mul.z*this->point3.z;
        return d;
    }
    void draw()
    {
        glColor3f(color[0],color[1],color[2]);
        glBegin(GL_TRIANGLES);
        {
            glVertex3f(point1.x, point1.y,point1.z);
            glVertex3f(point2.x,point2.y,point2.z);
            glVertex3f(point3.x,point3.y,point3.z);
        }
        glEnd();
    }
    void print()
    {
        this->point1.printPoint();
        this->point2.printPoint();
        this->point3.printPoint();
    }
    point getNormal(point p,Ray *r)
    {

        point tri= point(a, b, c);
        tri=tri.Normalize(tri);
        double check=r->direction.Dot_Multiplication(tri);
        if(check>0)
        {
            tri=tri*-1.0;
        }
        return tri;
    }
    double getintersection(Ray *ray)
    {
        const float EPSILON = 0.0000001;
        point edge1=this->point2-this->point1;
        point edge2=this->point3-this->point1;
        point temp=CrossProduct(ray->direction,edge2);
        double h=edge1.Dot_Multiplication(temp);
        if(h>-EPSILON && h<EPSILON)
        {
            return -1.0;
        }
        double fo=1.0/h;
        point s=ray->initial-this->point1;
        double u=fo*s.Dot_Multiplication(temp);
        if (u < 0.0 || u > 1.0)
        {
            return -1.0;
        }
        point q =CrossProduct(s,edge1);
        double v=fo*ray->direction.Dot_Multiplication(q);
        if (v < 0.0 || u + v > 1.0)
        {

            return -1.0;
        }
        double t = fo *edge2.Dot_Multiplication(q);
        //  cout<<"t......."<<t<<endl;
        if (t > EPSILON) // ray intersection
        {
            return t;
        }
        else // This means that there is a line intersection but not a ray intersection.
        {
            return -1;
        }


    }
    double intersection(Ray *r, double *new_color, int level)
    {

        double t=getintersection(r);
        if(t<=0 || isnan(t))
        {
            return -1;
        }
        if(level==0)
        {
            return t;
        }
        for(int i=0; i<3; i++)
        {
            new_color[i]=this->color[i]*this->coEfficient[0];
             if(new_color[i]<0.0)
            {
                new_color[i]=0.0;
            }
            else if(new_color[i]>1.0)
            {
                new_color[i]=1.0;
            }


        }

        point temp=r->direction*t;
        point intersectionPoint=r->initial+temp;
        point normal=getNormal(intersectionPoint,r);
        point reflection=getReflection(r,normal);
        double min_t=10000000;
        for(int i=0; i<light.size(); i++)
        {
            point dir=light[i]->ref_point-intersectionPoint;
            point norm_dir=dir.Normalize(dir);
            point initial=intersectionPoint+norm_dir*0.000001;
            Ray ro(initial,dir);
            bool flag=false;
            double t_min = 10000000;
            for(int j=0; j<object.size(); j++)
            {
                double t_temp=object[j]->getintersection(&ro);
                if(t_temp>0 && t_temp<t_min)
                {
                    t_min=t_temp;
                    flag=true;
                    break;

                }
            }
            if(!flag)
            {
                point d=ro.direction*-1.0;
                double dotn=d.Dot_Multiplication(normal)*2.0;
                point R=d-normal*dotn;
                R=R.Normalize(R);
                point V= r->direction *-1.0;
                V=V.Normalize(V);
                double lambert=max(0.0,ro.direction.Dot_Multiplication(normal));
                double phong = max(0.0,R.Dot_Multiplication(V));
                phong=pow(phong,shine);
                for (int k=0; k<3; k++)
                {
                    new_color[k] =new_color[k]+ light[i]->colors[k]* lambert * coEfficient[1] *this->color[k];
                    new_color[k]=new_color[k]+light[i]->colors[k] * phong* coEfficient[2] *this->color[k];
                    if(new_color[k]<0.0)
                    {
                        new_color[k]=0.0;
                    }
                    else if(new_color[k]>1.0)
                    {
                        new_color[k]=1.0;
                    }
                    // cout<<"new_color_now.............."<<new_color[k]<<endl;
                }

            }

        }
        if(level<recursion_level)
        {
             point initial=intersectionPoint+reflection*0.000001;
            Ray ref_ray(initial,reflection);
            int nearest=-1;
            double *ref_color=new double[3];

            for(int k=0; k<object.size(); k++)
            {
                double to = object[k]->getintersection(&ref_ray);
                if(to>0 && to<min_t )
                {
                    min_t =to;
                    nearest=k;
                }
            }
            if(nearest!=-1)
            {
                min_t = object[nearest]->intersection(&ref_ray,ref_color,level+1);
                for (int k=0; k<3; k++)
                {
                    new_color[k]+=ref_color[k]*coEfficient[3];
                     if(new_color[k]<0.0)
                    {
                        new_color[k]=0.0;
                    }
                    else if(new_color[k]>1.0)
                    {
                        new_color[k]=1.0;
                    }
                }
            }

        }
        return  t;
    }
};
class General:public Object
{
public:
    double A,B,C,D,E,F,G,H,I,J;
    point reference_point;
    double height,width,length;
    General(double A,double B,double C,double D,double E,double F,double G,double H,double I,double J,point ref_point,double length,double width,double height)
    {
        this->A=A;
        this->B=B;
        this->C=C;
        this->D=D;
        this->E=E;
        this->F=F;
        this->G=G;
        this->H=H;
        this->I=I;
        this->J=J;
        this->reference_point=ref_point;
        this->height=height;
        this->width=width;
        this->length=length;

    }
    void draw()
    {
    }
    void print()
    {
    }
    point getNormal(point p,Ray *r)
    {
        double a,b,c;
        a = 2 * this->A *p.x + this->D * p.y + this->F * p.z + G;
        b= 2 * this->B * p.y +this-> D * p.x +this-> E * p.z + H;
        c= 2 * this->C * p.z + this->E * p.y + this->F * p.x + I;

        point res;
        res.x=a;
        res.y=b;
        res.z=c;
        res=res.Normalize(res);
        double check=r->direction.Dot_Multiplication(res);
        if(check>0)
        {
            res=res*-1.0;
        }
        return res;

    }
    double getD(point p)
    {

    }
    double getintersection(Ray *ray)
    {
        double a=A*pow(ray->direction.x,2)+B*pow(ray->direction.y,2)+C*pow(ray->direction.z,2);
        double b=2*(A*ray->initial.x*ray->direction.x+B*ray->initial.y*ray->direction.y+C*ray->initial.z*ray->direction.z);
        double c=A*pow(ray->initial.x,2)+B*pow(ray->initial.y,2)+C*pow(ray->initial.z,2);
        a=a+(D*ray->direction.x*ray->direction.y+E*ray->direction.y*ray->direction.z+F*ray->direction.x*ray->direction.z);
        b=b+(D*(ray->initial.x*ray->direction.y+ray->direction.x*ray->initial.y)+E*(ray->initial.y*ray->direction.z+ray->direction.y*ray->initial.z)+F*((ray->initial.x*ray->direction.z+ray->direction.z*ray->initial.x)));
        c=c+(D*ray->initial.x*ray->initial.y+E*ray->initial.y*ray->initial.z+F*ray->initial.z*ray->initial.x);
        b=b+(G*ray->direction.x+H*ray->direction.y+I*ray->direction.z);
        c=c+(G*ray->initial.x+H*ray->initial.y+I*ray->initial.z+J);
        double dis=b*b-4*a*c;
        if(dis<0.0)
        {
            return -1;
        }
        else
        {
            double t1=(-b + sqrt(dis))/(2.0 * a);
            double t2 = (-b - sqrt(dis))/(2.0 * a);

            //  cout<<"t1.........."<<t1<<"t2............."<<t2<<endl;
            point intesec1=ray->initial+ray->direction*t1;
            point intesec2=ray->initial+ray->direction*t2;
            double min_x,max_x,min_y,max_y,min_z,max_z;
            min_x=this->reference_point.x;
            max_x=min_x+this->length;
            min_y=this->reference_point.y;
            max_y=min_y+this->width;
            min_z=this->reference_point.z;
            max_z=min_z+this->height;
            double a1,b1,c1,a2,b2,c2;
            a1=intesec1.x;
            b1=intesec1.y;
            c1=intesec1.z;
            a2=intesec2.x;
            b2=intesec2.y;
            c2=intesec2.z;
            bool f1=false;
            bool f2=false;
            if((this->length>0 && (a1<min_x || a1>max_x))||(this->width>0 &&(b1<min_y || b1>max_y))||(this->height>0 &&(c1<min_z|| c1>max_z)))
            {
                //cout<<"a....."<<a<<"min_x,max_x..."<<min_x<<max_x<<"b...."<<b<<"min_b,max_b..."<<min_y<<max_y<<"c..."<<c<<"minc maxc...."<<min_z<<max_z<<endl;
                f1=true;
            }
            if((this->length>0 && (a2<min_x || a2>max_x))||(this->width>0 &&(b2<min_y || b2>max_y))||(this->height>0 &&(c2<min_z || c2>max_z)))
            {
                f2=true;
            }
            if(f1==true && f2==true)
            {
                return -1;
            }
            else if(f1==true && f2==false)
            {
                //  cout<<"t2---"<<t2<<endl;

                return t2;
            }
            else if(f1==false && f2==true)
            {
                // cout<<"t2---"<<t1<<endl;
                return t1;
            }
            else
            {
                double res=min(t1,t2);
                return res;
            }


        }

    }
    double intersection(Ray *r, double *new_color, int level)
    {

        double t=getintersection(r);
        if(t<0)
        {
            return -1;
        }
        if(level==0)
        {
            return t;
        }
        for(int i=0; i<3; i++)
        {
            new_color[i]=this->color[i]*this->coEfficient[0];

        }
        //   cout<<"ashche............"<<endl;
        point temp=r->direction*t;
        point intersectionPoint=r->initial+temp;
        point normal=getNormal(intersectionPoint,r);
        point reflection=getReflection(r,normal);
        double min_t=10000000;
        for(int i=0; i<light.size(); i++)
        {
            point dir=light[i]->ref_point-intersectionPoint;
            point norm_dir=dir.Normalize(dir);
            point initial=intersectionPoint+norm_dir*0.000001;
            Ray ro(initial,norm_dir);
            bool flag=false;
            for(int j=0; j<object.size(); j++)
            {
                double t_temp=object[j]->getintersection(&ro);
                if(t_temp>0)
                {
                    flag=true;
                    break;

                }
            }
            if(!flag)
            {
                //lambert er jonne
                point d=ro.direction*-1.0;
                double dotn=d.Dot_Multiplication(normal)*2.0;
                point R=d-normal*dotn;
                R=R.Normalize(R);
                point V= r->direction *-1.0;
                V=V.Normalize(V);
                double lambert=max(0.0,ro.direction.Dot_Multiplication(normal));
                double phong = max(0.0,R.Dot_Multiplication(V));
                phong=pow(phong,shine);
                for (int k=0; k<3; k++)
                {
                    new_color[k] =new_color[k]+ light[i]->colors[k]* lambert * coEfficient[1] *this->color[k];
                    new_color[k]=new_color[k]+light[i]->colors[k] * phong * coEfficient[2] *this->color[k];
                    if(new_color[k]<0.0)
                    {
                        new_color[k]=0.0;
                    }
                    else if(new_color[k]>1.0)
                    {
                        new_color[k]=1.0;
                    }
                    // cout<<"new_color_now.............."<<new_color[k]<<endl;
                }

            }

        }
        if(level<recursion_level)
        {
            point initial=intersectionPoint+reflection*0.000001;
            Ray ref_ray(initial,reflection);
            int nearest=-1;

            double *ref_color=new double[3];
            for(int k=0; k<object.size(); k++)
            {

                double to = object[k]->getintersection(&ref_ray);
                if(to>0 && to<min_t )
                {
                    min_t =to;
                    nearest=k;
                }
            }
            if(nearest!=-1)
            {
                min_t = object[nearest]->intersection(&ref_ray,ref_color,level+1);
                for (int k=0; k<3; k++)
                {

                    new_color[k]+=ref_color[k]*coEfficient[3];
                     if(new_color[k]<0.0)
                    {
                        new_color[k]=0.0;
                    }
                    else if(new_color[k]>1.0)
                    {
                        new_color[k]=1.0;
                    }
                }
            }

        }
        return  t;



    }


};
class checkerBoard:public Object
{
public:
    double boardwidth,tileWidth,length;
    point ref_point;
    checkerBoard(double w,double t)
    {
        this->boardwidth=w;
        this->tileWidth=t;
        ref_point.x=-w/2;
        ref_point.y=-w/2;
        ref_point.z=0;

    }
    void draw()
    {

        int tiles_no=this->boardwidth/this->tileWidth;
        for(int i=0; i< tiles_no; i++)
        {
            for(int j=0; j<tiles_no; j++)
            {
                if ((i + j) % 2)
                {

                    glColor3f(0.0, 0.0, 0.0);
                }
                else
                {

                    glColor3f(1.0, 1.0, 1.0);
                }

                glBegin(GL_QUADS);
                {
                    glVertex3f(this->ref_point.x + this->tileWidth * i, this->ref_point.y + this->tileWidth * j, this->ref_point.z);
                    glVertex3f(this->ref_point.x + this->tileWidth * (i+1), this->ref_point.y + this->tileWidth * j, this->ref_point.z);
                    glVertex3f(this->ref_point.x + this->tileWidth *(i+1), this->ref_point.y + this->tileWidth * (j+1), this->ref_point.z);
                    glVertex3f(this->ref_point.x + this->tileWidth * i, this->ref_point.y + this->tileWidth * (j+1), this->ref_point.z);

                }
                glEnd();
            }

        }
    }
    void print()
    {

    }
    point getNormal(point p,Ray *r)
    {
        point ps;
        ps.x=0;
        ps.y=0;
        ps.z=1;
        double check=r->direction.Dot_Multiplication(ps);
        if(check>0)
        {
            ps=ps*-1.0;
        }
        return ps;
    }
    double getD(point p)
    {

    }
    double getintersection(Ray *ray)
    {
        double t;
        point start=ray->initial-this->ref_point;
        point normal;
        normal.x=0;
        normal.y=0;
        normal.z=1;
        double dn=normal.Dot_Multiplication(ray->direction);
        if(dn!=0)
        {
            t=-start.Dot_Multiplication(normal)/dn;
            return t;

        }
        else
        {
            t=-1.0;
        }
        return t;

    }
    double intersection(Ray *r, double *new_color, int level)
    {
        //   cout<<"z er val koto .............."<<this->ref_point.z<<endl;
        double t=this->getintersection(r);

        if(t<=0.0)
        {

            return -1.0;
        }
        if (level == 0)
        {
            return t;
        }

        point intersectionP=r->initial+r->direction*t;
        point new_p=intersectionP-this->ref_point;

        int new_I = new_p.x;
        int new_J = new_p.y;
        if(new_p.x >= this->boardwidth || new_p.y >= this->boardwidth || new_p.x < 0 || new_p.y < 0)
        {
            return -1.0;
        }

        for(int j = 0; j < abs(this->ref_point.y) * 2; j+=this->tileWidth)
        {
            for(int i = 0; i < abs(this->ref_point.x) * 2; i+=this->tileWidth)
            {

                if((new_I>=i && new_I<i+this->tileWidth) && (new_J >= j && new_J < j+this->tileWidth))
                {
                    if((i/20 + j/20) % 2 != 0)
                    {
                        new_color[0] = 0;
                        new_color[1] = 0;
                        new_color[2] = 0;


                    }
                    else
                    {
                        new_color[0] = 1;
                        new_color[1] = 1;
                        new_color[2] = 1;
                    }

                    break;
                }

            }
        }
        double *temp_color;
        temp_color=new double[3];
        for(int i=0;i<3;i++){
            temp_color[i]=new_color[i];
        }
        for(int i=0; i<3; i++)
        {
            new_color[i]=new_color[i]*0.6;
            //cout<<"new_color_before.............."<<new_color[i]<<endl;

        }

        point temp=r->direction*t;
        point intersectionPoint=r->initial+temp;
        point normal=getNormal(intersectionPoint,r);
        point reflection=getReflection(r,normal);
        double min_t=10000000;
        for(int i=0; i<light.size(); i++)
        {
            point dir=light[i]->ref_point-intersectionPoint;
            point norm_dir=dir.Normalize(dir);
            point initial=intersectionPoint+norm_dir*0.000001;
            Ray ro(initial,dir);
            bool flag=false;
            double t_min = 10000000;
            for(int j=0; j<object.size(); j++)
            {
                double t_temp=object[j]->getintersection(&ro);
                if(t_temp>0 && t_temp<t_min)
                {
                    t_min=t_temp;
                    flag=true;
                    break;

                }
            }
            if(!flag)
            {
                point d=ro.direction*-1.0;
                double dotn=d.Dot_Multiplication(normal)*2.0;
                point R=d-normal*dotn;
                R=R.Normalize(R);
                point V= r->direction *-1.0;
                V=V.Normalize(V);
                double lambert=max(0.0,ro.direction.Dot_Multiplication(normal));
                double phong = max(0.0,R.Dot_Multiplication(V));
                phong=pow(phong,1.0);
                //  cout<<"after power...."<<phong<<endl;


                for (int k=0; k<3; k++)
                {

                    new_color[k] =new_color[k]+ light[i]->colors[k]* lambert * 0.6 *temp_color[k];
                    new_color[k]=new_color[k]+light[i]->colors[k] *phong * 0.6 *temp_color[k];

                    if(new_color[k]<0.0)
                    {
                        new_color[k]=0.0;
                    }
                    else if(new_color[k]>1.0)
                    {
                        new_color[k]=1.0;
                    }
                    // cout<<"new_color_now.............."<<new_color[k]<<endl;
                }

            }

        }
        if(level<recursion_level)
        {
            point initial=intersectionPoint+reflection*0.000001;
            Ray ref_ray(initial,reflection);
            int nearest=-1;

            double *ref_color=new double[3];
            double min_t=1000000;
            for(int k=0; k<object.size(); k++)
            {


                double to = object[k]->getintersection(&ref_ray);
                if(to>0 && to<min_t )
                {
                    min_t =to;
                    nearest=k;
                }


            }
            if(nearest!=-1)
            {
                min_t = object[nearest]->intersection(&ref_ray,ref_color,level+1);
                for (int k=0; k<3; k++)
                {
                    new_color[k]+=ref_color[k]*0.6;
                     if(new_color[k]<0.0)
                    {
                        new_color[k]=0.0;
                    }
                    else if(new_color[k]>1.0)
                    {
                        new_color[k]=1.0;
                    }
                }
            }


        }
        return  t;


    }

};

#endif
