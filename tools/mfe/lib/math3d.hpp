#ifndef MATH3D_H
#define MATH3D_H

#include <cmath>

namespace Math3D {
    
    #define SINE_FUNCTION std::sin
    #define COSINE_FUNCTION std::cos
    #define SQRT_FUNCTION std::sqrt
    
    template<class T> class Vector4;
    template<class T> class Matrix4;
    
    template<class T> class Vector4 {
        public:
        Vector4(void) {}
        Vector4(T a, T b, T c, T d = 1)
        { _v[0]=a; _v[1]=b; _v[2]=c; _v[3]=d; }
        
        // The int parameter is the number of elements to copy from initArray (3 or 4)
        explicit Vector4(T* initArray, int arraySize = 3)
        { for (int i = 0;i<arraySize;++i) _v[i] = initArray[i]; if (arraySize == 3) _v[3] = 1; }
        
        // [] is to read, () is to write (const correctness)
        const T& operator[] (int i) const { return _v[i]; }
        T& operator() (int i) { return _v[i]; }
        
        inline bool operator== (const Vector4<T>&);
        inline bool operator< (const Vector4<T>&);
        
        inline Vector4<T>   operator+   (const Vector4<T>&);
        inline Vector4<T>&  operator+=  (const Vector4<T>&);
        inline Vector4<T>   operator-   (const Vector4<T>&);
        inline Vector4<T>&  operator-=  (const Vector4<T>&);
        inline Vector4<T>   operator*   (T);
        inline Vector4<T>&  operator*=  (T);
        inline Vector4<T>   operator/   (T);
        inline Vector4<T>&  operator/=  (T);

        inline T LengthSquared3();
        inline T Length3();
        inline Vector4<T>& Normalize3();
        inline T DotProduct3 (const Vector4<T>& _v);
        inline T CrossProduct (const Vector4<T>&, const Vector4<T>&);

        // Provides access to the underlying array; useful for passing this class off to C APIs
        const T* readArray(void) { return _v; }
        T* getArray(void) { return _v; }

        private:
        T _v[4];
    };

    template<class T> class Matrix4 {
        public:
        Matrix4 (void) {}
        
        
        // [] is to read, () is to write (const correctness)
        // m[x][y] or m(x)[y] is the correct form
        const T* operator[] (int i) const { return &m[i<<2]; }
        T* operator() (int i) { return &m[i<<2]; }
        
        // Low-level access to the array.
        const T* readArray (void) { return m; }
        T* getArray(void) { return m; }
        
        // Construct various matrices; REPLACES CURRENT CONTENTS OF THE MATRIX!
        // Written this way to work in-place and hence be somewhat more efficient
        void Identity (void) { for (int i=0;i<16;++i) m[i] = 0; m[0] = 1; m[5] = 1; m[10] = 1; m[15] = 1; }
        inline Matrix4& Rotation (T angle, Vector4<T> axis);
        inline Matrix4& Translation(const Vector4<T>& translation);
        inline Matrix4& Scale (T x, T y, T z);
        inline Matrix4& BasisChange (const Vector4<T>& v, const Vector4<T>& n);
        inline Matrix4& BasisChange (const Vector4<T>& u, const Vector4<T>& v, const Vector4<T>& n);
        inline Matrix4& ProjectionMatrix (bool perspective, T l, T r, T t, T b, T n, T f);
        
        private:
        T m[16];
    };
    
    template<class T> class Line4 {  /* use one point and the direction to record a line */
        public:
        Line4(void) {}
        Line4(Vector4<T> a, Vector4<T> b) {v0 = a; n = (b-a).Normalize3();}

        private:
        Vector4<T> v0;
        Vector4<T> n;
    };

    template<class T> class Line4_PP {  /* use two points to record a line */
        public:
        Line4_PP(void) {}
        Line4_PP(Vector4<T> a, Vector4<T> b) {v0 = a; v1 = b;}
        
        private:
        Vector4<T> v0;
        Vector4<T> v1;
    };
    
    /* Vector operators */
    template<class T> Vector4<T> Vector4<T>::operator+ (const Vector4<T>& v)
    {Vector4<T> v_new; for (int i=0;i<4;++i) v_new(i) = _v[i] + v[i]; return v_new; }
    
    template<class T> Vector4<T> Vector4<T>::operator- (const Vector4<T>& v)
    {Vector4<T> v_new; for (int i=0;i<4;++i) v_new(i) = _v[i] - v[i]; return v_new; }
    
    template<class T> Vector4<T>& Vector4<T>::operator+= (const Vector4<T>& v)
    { for (int i=0;i<4;++i) _v[i] += v[i]; return *this; }
    
    template<class T> Vector4<T>& Vector4<T>::operator-= (const Vector4<T>& v)
    { for (int i=0;i<4;++i) _v[i] -= v[i]; return *this; }
    
    template<class T> Vector4<T> Vector4<T>::operator* (const T c)
    {Vector4<T> v_new; for (int i=0;i<4;++i) v_new(i) = _v[i] * c; return v_new; }
    
    template<class T> Vector4<T>& Vector4<T>::operator*= (T c)
    { for (int i=0;i<4;++i) _v[i] *= c; return *this; }
    
    template<class T> Vector4<T>& Vector4<T>::operator/= (T c)
    { for (int i=0;i<4;++i) _v[i] /= c; return *this; }
    
    template<class T> T Vector4<T>::LengthSquared3()
    { return this->DotProduct3(*this); }

    template<class T> T Vector4<T>::Length3()
    { return SQRT_FUNCTION(this->LengthSquared3()); }
    
    template<class T> Vector4<T>& Vector4<T>::Normalize3()
    {*this /= this->Length3(); return *this;}

    template<class T> T Vector4<T>::DotProduct3(const Vector4<T>& v)
    { return _v[0]*v[0] + _v[1]*v[1] + _v[2]*v[2]; }
    
    template<class T> Vector4<T> CrossProduct (const Vector4<T>& v1, const Vector4<T>& v2) {
        return Vector4<T>( 	 v1[1] * v2[2] - v1[2] * v2[1]
        ,v2[0] * v1[2] - v2[2] * v1[0]
        ,v1[0] * v2[1] - v1[1] * v2[0]
        ,1);
    }

}

#endif

