//
//  This file is part of EPOS4
//  Copyright (C) 2022 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License version 3 or later
//  (See COPYING file for the text of the licence)
//


//###########################################################################

template <typename T>
class Mudiar
{
private:
	T *mudiar ;
        int m1, m2, m3, m4, m5;
public:
	Mudiar(int n1);
	Mudiar(int n1, int n2);
	Mudiar(int n1, int n2, int n3);
	Mudiar(int n1, int n2, int n3, int n4);
	Mudiar(int n1, int n2, int n3, int n4, int n5);
	~Mudiar(){ delete [] mudiar ; };
        void set(int i1, T val);
        void get(int i1, T &val);
        void set2(T v0,T v1);
        void get2(T &v0,T &v1);
        void set3(T v0,T v1,T v2);
        void get3(T &v0,T &v1,T &v2);
        void set4(T v0,T v1,T v2,T v3);
        void get4(T &v0,T &v1,T &v2,T &v3);
        void set5(T v0,T v1,T v2,T v3,T v4);
        void get5(T &v0,T &v1,T &v2,T &v3,T &v4);
        void set6(T v0,T v1,T v2,T v3,T v4,T v5);
        void get6(T &v0,T &v1,T &v2,T &v3,T &v4,T &v5);
        void set7(T v0,T v1,T v2,T v3,T v4,T v5,T v6);
        void get7(T &v0,T &v1,T &v2,T &v3,T &v4,T &v5,T &v6);
        void acc(int i1, T val);
        void zero();
        void set(int i1, int i2, T val);
        void get(int i1, int i2, T &val);
        void acc(int i1, int i2, T val);
        void increment(int i1, int i2);
        void set(int i1, int i2, int i3, T val);
        void get(int i1, int i2, int i3, T &val);
        void acc(int i1, int i2, int i3, T val);
        void set(int i1, int i2, int i3, int i4, T val);
        void get(int i1, int i2, int i3, int i4, T &val);
        void acc(int i1, int i2, int i3, int i4, T val);
        void set(int i1, int i2, int i3, int i4, int i5, T val);
        void get(int i1, int i2, int i3, int i4, int i5, T &val);
        void increment(int i1, int i2, int i3, int i4, int i5);
};
//------------------------------------1d-----------------------------
template<typename T>
Mudiar<T>::Mudiar(int n1){
	mudiar = new T [n1] ;
	for(int i=0; i<n1; i++){ mudiar[i]=0; }
        m1=n1; m2=1; m3=1; m4=1; m5=1;
}
template<typename T>
void Mudiar<T>::set(int i1, T val)        {  mudiar[i1] = val;  }
template<typename T>
void Mudiar<T>::get(int i1, T &_val)      { _val = mudiar[i1]; }
template<typename T>
void Mudiar<T>::set2(T v0,T v1)      {  mudiar[0] = v0; mudiar[1] = v1;  }
template<typename T>
void Mudiar<T>::get2(T &_v0,T &_v1)      { _v0 = mudiar[0]; _v1 = mudiar[1]; }
template<typename T>
void Mudiar<T>::set3(T v0,T v1,T v2)      {  mudiar[0] = v0; mudiar[1] = v1; mudiar[2] = v2;  }
template<typename T>
void Mudiar<T>::get3(T &_v0,T &_v1,T &_v2)      { _v0 = mudiar[0]; _v1 = mudiar[1]; _v2 = mudiar[2]; }
template<typename T>
void Mudiar<T>::set4(T v0,T v1,T v2,T v3)      {mudiar[0]=v0;mudiar[1]=v1;mudiar[2]=v2;mudiar[3]=v3;}
template<typename T>
void Mudiar<T>::get4(T &_v0,T &_v1,T &_v2,T &_v3)   {_v0=mudiar[0];_v1=mudiar[1];_v2=mudiar[2];_v3=mudiar[3];}
template<typename T>
void Mudiar<T>::set5(T v0,T v1,T v2,T v3,T v4)      {mudiar[0]=v0;mudiar[1]=v1;mudiar[2]=v2;mudiar[3]=v3;mudiar[4]=v4;}
template<typename T>
void Mudiar<T>::get5(T &_v0,T &_v1,T &_v2,T &_v3,T &_v4)   {_v0=mudiar[0];_v1=mudiar[1];_v2=mudiar[2];_v3=mudiar[3];_v4=mudiar[4];}
template<typename T>
void Mudiar<T>::set6(T v0,T v1,T v2,T v3,T v4,T v5)      {mudiar[0]=v0;mudiar[1]=v1;mudiar[2]=v2;mudiar[3]=v3;mudiar[4]=v4;mudiar[5]=v5;}
template<typename T>
void Mudiar<T>::get6(T &_v0,T &_v1,T &_v2,T &_v3,T &_v4,T &_v5)   {_v0=mudiar[0];_v1=mudiar[1];_v2=mudiar[2];_v3=mudiar[3];_v4=mudiar[4];_v5=mudiar[5];}
template<typename T>
void Mudiar<T>::set7(T v0,T v1,T v2,T v3,T v4,T v5,T v6)      {mudiar[0]=v0;mudiar[1]=v1;mudiar[2]=v2;mudiar[3]=v3;mudiar[4]=v4;mudiar[5]=v5;mudiar[6]=v6;}
template<typename T>
void Mudiar<T>::get7(T &_v0,T &_v1,T &_v2,T &_v3,T &_v4,T &_v5,T &_v6)   {_v0=mudiar[0];_v1=mudiar[1];_v2=mudiar[2];_v3=mudiar[3];_v4=mudiar[4];_v5=mudiar[5];_v6=mudiar[6];}
template<typename T>
void Mudiar<T>::acc(int i1, T val)        {  mudiar[i1] = mudiar[i1-1] + val; }
template<typename T>
void Mudiar<T>::zero()                { for(int i=0; i<m1; i++){ mudiar[i]=0; } }
//------------------------------------2d-----------------------------
template<typename T>
Mudiar<T>::Mudiar(int n1, int n2){
	mudiar = new T [n1*n2] ;
	for(int i=0; i<n1*n2; i++){mudiar[i]=0;}
        m1=n1; m2=n2; m3=1; m4=1; m5=1;
}
template<typename T>
void Mudiar<T>::set(int i1, int i2, T val)  { mudiar[i1+m1*i2] = val; }
template<typename T>
void Mudiar<T>::get(int i1, int i2, T &_val){ _val = mudiar[i1+m1*i2]; }
template<typename T>
void Mudiar<T>::acc(int i1, int i2, T val)  { mudiar[i1+m1*i2] = mudiar[i1+m1*(i2-1)] + val; }
template<typename T>
void Mudiar<T>::increment(int i1, int i2)   { mudiar[i1+m1*i2] = mudiar[i1+m1*i2] + 1 ; }
//------------------------------------3d-----------------------------
template<typename T>
Mudiar<T>::Mudiar(int n1, int n2, int n3){
	mudiar = new T [n1*n2*n3] ;
	for(int i=0; i<n1*n2*n3; i++){mudiar[i]=0;}
        m1=n1; m2=n2; m3=n3; m4=1; m5=1;
}
template<typename T>
void Mudiar<T>::set(int i1, int i2, int i3, T val)  { mudiar[i1+m1*i2+m1*m2*i3] = val; }
template<typename T>
void Mudiar<T>::get(int i1, int i2, int i3, T &_val){ _val = mudiar[i1+m1*i2+m1*m2*i3]; }
template<typename T>
void Mudiar<T>::acc(int i1, int i2, int i3, T val)  { mudiar[i1+m1*i2+m1*m2*i3] = mudiar[i1+m1*i2+m1*m2*(i3-1)] + val; }
//------------------------------------4d-----------------------------
template<typename T>
Mudiar<T>::Mudiar(int n1, int n2, int n3, int n4){
	mudiar = new T [n1*n2*n3*n4] ;
	for(int i=0; i<n1*n2*n3*n4; i++){mudiar[i]=0;}
        m1=n1; m2=n2; m3=n3; m4=n4; m5=1;
}
template<typename T>
void Mudiar<T>::set(int i1, int i2, int i3, int i4, T val)  { mudiar[i1+m1*i2+m1*m2*i3+m1*m2*m3*i4] = val; }
template<typename T>
void Mudiar<T>::get(int i1, int i2, int i3, int i4, T &_val){ _val = mudiar[i1+m1*i2+m1*m2*i3+m1*m2*m3*i4]; }
template<typename T>
void Mudiar<T>::acc(int i1, int i2, int i3, int i4, T val)  { mudiar[i1+m1*i2+m1*m2*i3+m1*m2*m3*i4] = mudiar[i1+m1*i2+m1*m2*i3+m1*m2*m3*(i4-1)] + val; }
//------------------------------------5d-----------------------------
template<typename T>
Mudiar<T>::Mudiar(int n1, int n2, int n3, int n4, int n5){
	mudiar = new T [n1*n2*n3*n4*n5] ;
	for(int i=0; i<n1*n2*n3*n4*n5; i++){mudiar[i]=0;}
        m1=n1; m2=n2; m3=n3; m4=n4; m5=n5;
}
template<typename T>
void Mudiar<T>::set(int i1, int i2, int i3, int i4, int i5, T val)  { mudiar[i1+m1*i2+m1*m2*i3+m1*m2*m3*i4+m1*m2*m3*m4*i5] = val; }
template<typename T>
void Mudiar<T>::get(int i1, int i2, int i3, int i4, int i5, T &_val){ _val = mudiar[i1+m1*i2+m1*m2*i3+m1*m2*m3*i4+m1*m2*m3*m4*i5]; }
template<typename T>
void Mudiar<T>::increment(int i1, int i2, int i3, int i4, int i5)   { mudiar[i1+m1*i2+m1*m2*i3+m1*m2*m3*i4+m1*m2*m3*m4*i5]++ ; }

//###########################################################################

class HoTab
{
private:
	float *velio, *bario ;
        int m1, m2, m3, m4, m5;
public:
	HoTab(int n1, int n2, int n3, int n4, int n5);
	~HoTab();
        void velioset(int i1, int i2, int i3, int i4, int i5, float val);
        void barioset(int i1, int i2, int i3, int i4, int i5, float val);
        void velioget(int i1, int i2, int i3, int i4, int i5, float &val);
        void barioget(int i1, int i2, int i3, int i4, int i5, float &val);
};

//###########################################################################

class CCCptl
{
private:
        int id, ist, ity, ior, jor, ifr1, ifr2 ;
        float tiv1, tiv2, p1, p2, p3,  p4, p5 ;
        float xor1, xor2, xor3, xor4, rad ;
        float des, dez, qsq, zpa1, zpa2, rin ;
        int ib1, ib2, ib3, ib4, iaa, its ; 
public:
        CCCptl(): id(0), ist(0), ity(0), ior(0), jor(0), ifr1(0), ifr2(0),
             tiv1(0), tiv2(0), p1(0), p2(0), p3(0),  p4(0), p5(0),
             xor1(0), xor2(0), xor3(0), xor4(0), rad(0),
             des(0), dez(0), qsq(0), zpa1(0), zpa2(0), rin(0),
             ib1(0), ib2(0), ib3(0), ib4(0), iaa(0), its(0) {} ;
        ~CCCptl() {}; // empty destructor, since nothing is allocated in heap (dynamic vars, etc)
        void dump(    int _id, int _ist, int _ity, int _ior, int _jor, int _ifr1, int _ifr2
                    , float _tiv1, float _tiv2, float _p1, float _p2, float _p3, float _p4, float _p5
                    , float _xor1, float _xor2, float _xor3, float _xor4, float _rad
                    , float _des, float _dez, float _qsq, float _zpa1, float _zpa2, float _rin
                    , int _ib1, int _ib2, int _ib3, int _ib4, int _iaa, int _its) 
                 {
                    id    = _id    ;
                    ist   = _ist   ;
                    ity   = _ity   ;
                    ior   = _ior   ;
                    jor   = _jor   ;
                    ifr1  = _ifr1  ;
                    ifr2  = _ifr2  ;
            
                    tiv1  = _tiv1  ;
                    tiv2  = _tiv2  ;
                    p1    = _p1    ;
                    p2    = _p2    ;
                    p3    = _p3    ;
                    p4    = _p4    ;
                    p5    = _p5    ;
            
                    xor1  = _xor1  ;
                    xor2  = _xor2  ;
                    xor3  = _xor3  ;
                    xor4  = _xor4  ;
                    rad   = _rad   ;
            
                    des   = _des   ;
                    dez   = _dez   ;
                    qsq   = _qsq   ;
                    zpa1  = _zpa1  ;
                    zpa2  = _zpa2  ;
                    rin   = _rin   ;
            
                    ib1   = _ib1   ;
                    ib2   = _ib2   ;
                    ib3   = _ib3   ;
                    ib4   = _ib4   ;
                    iaa   = _iaa   ;
                    its   = _its   ;
                 }

        void restore( int &_id, int &_ist, int &_ity, int &_ior, int &_jor, int &_ifr1, int &_ifr2
                    , float &_tiv1, float &_tiv2, float &_p1, float &_p2, float &_p3, float &_p4, float &_p5
                    , float &_xor1, float &_xor2, float &_xor3, float &_xor4, float &_rad
                    , float &_des, float &_dez, float &_qsq, float &_zpa1, float &_zpa2, float &_rin
                    , int &_ib1, int &_ib2, int &_ib3, int &_ib4, int &_iaa, int &_its) const
                 {
                    _id    = id    ;
                    _ist   = ist   ;
                    _ity   = ity   ;
                    _ior   = ior   ;
                    _jor   = jor   ;
                    _ifr1  = ifr1  ;
                    _ifr2  = ifr2  ;
            
                    _tiv1  = tiv1  ;
                    _tiv2  = tiv2  ;
                    _p1    = p1    ;
                    _p2    = p2    ;
                    _p3    = p3    ;
                    _p4    = p4    ;
                    _p5    = p5    ;
            
                    _xor1  = xor1  ;
                    _xor2  = xor2  ;
                    _xor3  = xor3  ;
                    _xor4  = xor4  ;
                    _rad   = rad   ;
            
                    _des   = des   ;
                    _dez   = dez   ;
                    _qsq   = qsq   ;
                    _zpa1  = zpa1  ;
                    _zpa2  = zpa2  ;
                    _rin   = rin   ;
                    
                    _ib1   = ib1   ;
                    _ib2   = ib2   ;
                    _ib3   = ib3   ;
                    _ib4   = ib4   ;
                    _iaa   = iaa   ;
                    _its   = its   ;
                 }

        inline void  setID (int value) { id = value ; }
        inline void setIST (int value) { ist = value ; }
        inline void setITY (int value) { ity = value ; }
        inline void setIOR (int value) { ior = value ; } 
        inline void setJOR (int value) { jor = value ; } 
        inline void setIFR (int _i1, int _i2) {ifr1=_i1;ifr2=_i2;} 
        inline void setDES (float value) { des = value ; }
        inline void setRAD (float value) { rad = value ; }
        inline void setRIN (float value) { rin = value ; }
        inline void setP(float _p1, float _p2, float _p3, float _p4, float _p5) {p1=_p1;p2=_p2;p3=_p3;p4=_p4;p5=_p5;}
        inline void setXOR(float _x1, float _x2, float _x3, float _x4) {xor1=_x1;xor2=_x2;xor3=_x3;xor4=_x4;}
        inline void setTIV(float _t1, float _t2) {tiv1=_t1;tiv2=_t2;}
        inline void setZPA(float _t1, float _t2) {zpa1=_t1;zpa2=_t2;}
        inline void setIB (int _i1, int _i2, int _i3, int _i4) {ib1=_i1;ib2=_i2;ib3=_i3;ib4=_i4;} 
        inline void set2IB (int j, int _ival) {if(j==1){ib1=_ival;}if(j==2){ib2=_ival;}if(j==3){ib3=_ival;}if(j==4){ib4=_ival;}} 

        inline int  getID () const {return id ; }
        inline int getIST () const {return ist ; }
        inline int getITY () const {return ity ; }  
        inline int getIOR () const {return ior ; } 
        inline int getJOR () const {return jor ; } 
        inline void getIFR (int &_i1, int &_i2) const {_i1=ifr1;_i2=ifr2;} 
        inline float getDES () const {return des ; }
        inline float getRAD () const {return rad ; }
        inline float getRIN () const {return rin ; }
        inline void getP(float &_p1, float &_p2, float &_p3, float &_p4, float &_p5) const {_p1=p1;_p2=p2;_p3=p3;_p4=p4;_p5=p5;}
        inline void getXOR(float &_x1, float &_x2, float &_x3, float &_x4) const {_x1=xor1;_x2=xor2;_x3=xor3;_x4=xor4;}
        inline void getTIV(float &_t1, float &_t2) const {_t1=tiv1;_t2=tiv2;}
        inline void getZPA(float &_t1, float &_t2) const {_t1=zpa1;_t2=zpa2;}
        inline void getIB (int &_i1, int &_i2, int &_i3, int &_i4) const {_i1=ib1;_i2=ib2;_i3=ib3;_i4=ib4;} 
        inline void get2IB (int j, int &_ival) const {if(j==1){_ival=ib1;}if(j==2){_ival=ib2;}if(j==3){_ival=ib3;}if(j==4){_ival=ib4;}} 
};

//###########################################################################

   /* OmTab removed after version 3238 */



