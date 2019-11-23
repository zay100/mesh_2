//======================================================================================================================
// Model implementation of containers for mesh data
//======================================================================================================================
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <stdarg.h>
using namespace std;

//======================================================================================================================
// Global definitions 
//======================================================================================================================

#ifdef MPI_PARALLEL
    #include <mpi.h>
#else
    #define MPI_Comm int
    #define MPI_COMM_WORLD 0    
#endif

#define MAX_UINT64             (0xEFFFFFFFFFFFFFFFull)
#define MAX_ALLOWED_BLOCK_SIZE (0xEFFFFFFFl) // 4026531839(3.8Gb)
#define MAX_SIGNED_INT         (0x7FFFFFFF ) // 2147483647 

#define ULL(X) ((unsigned long long)X) // use %llu to print 

#ifndef SAFE_MODE
#define SAFE_MODE
#endif

// Assert macro wrappers 
#ifdef SAFE_MODE
// to be placed in performance-critical places - enabled in safe-mode build only
#define SAFE_ASSERT(X,...) if(!(X)) exit(Crash((stringprintf("Assert (%s) failed: ", #X) + stringprintf(__VA_ARGS__)).c_str()));
#else
#define SAFE_ASSERT(X,...) // disabled
#endif
// always enabled - to be placed at upper level (where check doesn't affect performance)
#define ASSERT(X,...)      if(!(X)) exit(Crash((stringprintf("Assert (%s) failed: ", #X) + stringprintf(__VA_ARGS__)).c_str()));

#pragma warning (disable: 4127) // отключить ругань майкрософтовского компеллера с W4 на константы в условном операторе
#pragma warning (disable: 6993) // ворнинг аналайзера об игнорировании ОпенМП
#pragma warning (disable: 4068) // ворнинг аналайзера об игнорировании незнакомой прагмы ivdep
#pragma warning (disable: 6211) // ворнинг на каждый new, что надо ловить эксепшн. пока отключено. 

//======================================================================================================================
// Utility functions 
//======================================================================================================================

// -----------------------------------------------------------------------------------------------------------
// Global data
// -----------------------------------------------------------------------------------------------------------
int mpi_initialized = 0; // flag that MPI is initialized 
int MyID = 0; // process ID
int NumProc = 1; // number of processes 
int MASTER_ID = 0; // master process ID
MPI_Comm MCW = MPI_COMM_WORLD; // default communicator

#define LINE_SEPARATOR  "-------------------------------------------------------------------------------\n" 

// -----------------------------------------------------------------------------------------------------------
// Diagnostics 
// -----------------------------------------------------------------------------------------------------------

#define crash(...) exit(Crash( __VA_ARGS__)) // via exit define so static analyzer knows its an exit point
static int Crash(const char *fmt,...){ // termination of program due to error
    va_list ap;
    if(mpi_initialized) fprintf(stderr,"\nEpic fail: MyID = %d\n",MyID);
    else fprintf(stderr,"\nEpic fail: \n");

    va_start(ap,fmt);
    vfprintf(stderr,fmt,ap);
    va_end(ap);

    fprintf(stderr,"\n");
    fflush(stderr);

#ifdef MPI_PARALLEL    
    if(mpi_initialized) MPI_Abort(MPI_COMM_WORLD,-1);
#endif
    return 0;
}

// Internal interface to MPI barrier
static inline void barrier(){
#ifdef MPI_PARALLEL 
     if(MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) crash("Base lib barrier: MPI_Barrier failed! \n");
#endif
}

int printf0(const char *fmt,...){ // Write to stdout from Master process
    int r = 0;
    va_list ap;
    if(MyID==MASTER_ID){
        va_start(ap,fmt);  
        r=vfprintf(stdout,fmt,ap);  
        va_end(ap);
    }
    fflush(stdout);
    return(r);
}

// Debug synchronous printf - Write to stdout + flush + barrier
int pprintf(const char *fmt,...){ 
    int r = 0;
    fflush(stdout);
    barrier();
    for(int p=0; p<NumProc; p++){
        barrier();
        if(MyID != p) continue; 
        fprintf(stdout,"%3d: ",MyID);
        va_list ap;
        //stdout
        va_start(ap,fmt);  
        r=vfprintf(stdout,fmt,ap);
        va_end(ap);
        fflush(stdout);
    }
    fflush(stdout);
    barrier();
    return(r);
}

// Formatted output to std::string 
string stringprintf(const char* fmt, ...) {
    int bufsize = 128; // начальный размер буфера
    char *sbuf = new char [bufsize]; //буфер для пихания в чар-строку

    do{//пока успешно не принтанем 
        va_list ap;
        va_start(ap,fmt);
        int nstr = vsnprintf(sbuf,bufsize,fmt,ap);
        va_end(ap);

        if((nstr>=0)&&(nstr<bufsize)){//корректно принтанулось
            string s =  string(sbuf); //результирующий стринг
            delete[] sbuf; 
            return s;            
        }
        if(nstr >= bufsize) bufsize = nstr + 1;
        if(nstr < 0) bufsize *= 2;
        
        delete[] sbuf;
        sbuf = new char [bufsize];

    }while(bufsize < 1000000); //пока строка не разрослась до абсурдных размеров
    delete[] sbuf;
    crash("stringprintf failed! buffer size too big."); 
}

// -----------------------------------------------------------------------------------------------------------
// Memory allocation
// -----------------------------------------------------------------------------------------------------------
 
template <typename ValueType>
inline ValueType* GimmeMem(size_t N, const char* label = NULL){
    if(N==0) return NULL;

    if(sizeof(size_t)==sizeof(int)){// size check for 32-bit system
        if(N*sizeof(ValueType) > MAX_ALLOWED_BLOCK_SIZE) 
            crash("%s: size %llu exceeds maximal value %llu", 
                  label?label:"", ULL(N*sizeof(ValueType)), ULL(MAX_ALLOWED_BLOCK_SIZE));
    }

    ValueType* p=NULL;
    try{p = new ValueType[N];}
    catch (bad_alloc &ba){
        if(p!=NULL) delete[] p;
        crash("GimmeMem: memory fuckup! N = %llu, sizeof = %d, name %s, text: %s\n", 
              (unsigned long long) N, (int)sizeof(ValueType), label ? label:"noname",ba.what());
    }
    return p;
}

// deallocation wrapper
template <typename ValueType> 
inline void FreeMem(ValueType* __restrict &Array) { 
    if(!Array) crash("NULL DEALLOC!"); //на всякий случай. 
    { delete[] Array; Array=NULL; }
}
//======================================================================================================================


//======================================================================================================================
// Index types
//======================================================================================================================

typedef int tInt; // integer index value (to switch 32/64 bit if needed)

// Kinds of mesh elements
enum tIdxKind{ IDX_N=0/*nodes*/, IDX_E=1/*elements*/, IDX_S=2/*segments*/, IDX_F=3/*faces*/ /*, ...*/, IDX_C=4};

template<tIdxKind I> class tIdx { // Index of a given kind - wrapper around an integer value
private:    
    tInt i; // index value 
 //   inline tIdx(const tInt j){ i=j; } // forbids direct cast from ints 
    //  inline operator tInt() const {return i;}  // forbidden
public:
    inline tInt idx() const { return i;}
    inline tIdx(){ i=-1; } // an incorrect value by default 
    inline void Set(int j){i=j;}
    inline tIdx<I>& operator=(const tInt   j){   i=j;   return *this; }
    inline tIdx<I>& operator=(const size_t j){   i=(tInt)j;   return *this; }
    inline tIdx<I>& operator=(const tIdx<I>& j){ i=j.i; return *this; }

    inline tIdx<I>& operator++() { ++i; return *this; } 
    inline tIdx<I>& operator--() { --i; return *this; } 
    
    inline tIdx<I> operator++(int) { return tIdx<I>(i++); } 
    inline tIdx<I> operator--(int) { return tIdx<I>(i--); } 

    inline tIdx<I> operator+(tIdx<I> j) { return tIdx<I>(i+j.i); }
    inline tIdx<I> operator-(tIdx<I> j) { return tIdx<I>(i-j.i); }

    inline tIdx<I> operator+(tInt j) { return tIdx<I>(i+j); }
    inline tIdx<I> operator-(tInt j) { return tIdx<I>(i+j); }

    inline tIdx<I> operator*(tInt j) { return tIdx<I>(i*j); }
    inline tIdx<I> operator/(tInt j) { return tIdx<I>(i/j); }
    inline tIdx<I> operator%(tInt j) { return tIdx<I>(i%j); }

    inline bool operator< (tInt j) const { return i <  j; }
    inline bool operator> (tInt j) const { return i >  j; }
    inline bool operator<=(tInt j) const { return i <= j; }
    inline bool operator>=(tInt j) const { return i >= j; }
    inline bool operator==(tInt j) const { return i == j; }
    inline bool operator!=(tInt j) const { return i != j; }

    inline bool operator< (tIdx<I> j) const { return i <  j.i; }
    inline bool operator> (tIdx<I> j) const { return i >  j.i; }
    inline bool operator<=(tIdx<I> j) const { return i <= j.i; }
    inline bool operator>=(tIdx<I> j) const { return i >= j.i; }
    inline bool operator==(tIdx<I> j) const { return i == j.i; }
    inline bool operator!=(tIdx<I> j) const { return i != j.i; }

    inline bool operator!() { return !i; }

    template<tIdxKind TT> tIdx(const tIdx<TT> &J);
};
 
inline bool IsCC(){return false;} // Checks if the numerical scheme is cell-centered
template<>template<> inline tIdx<IDX_N>::tIdx(const tIdx<IDX_C>& j){SAFE_ASSERT(!IsCC(),"Wrong C2N"); i=j.idx();}
template<>template<> inline tIdx<IDX_C>::tIdx(const tIdx<IDX_N>& j){SAFE_ASSERT(!IsCC(),"Wrong N2C"); i=j.idx();}
template<>template<> inline tIdx<IDX_E>::tIdx(const tIdx<IDX_C>& j){SAFE_ASSERT( IsCC(),"Wrong C2E"); i=j.idx();}
template<>template<> inline tIdx<IDX_C>::tIdx(const tIdx<IDX_E>& j){SAFE_ASSERT( IsCC(),"Wrong E2C"); i=j.idx();}


//======================================================================================================================
// Basic containers 
// tBlock       - small fixed-size array (vector), used in block arrays
// tArray       - big allocatable array
// tBlockArray  - 2D array [big][small] - block vector with small blocks of equal size. :tArray.
// tVBlockArray - same as tBlockArray but with blocks of different size. :tArray. mainly used to store CSR topology
//======================================================================================================================
#define SAFE_MODE 

#ifdef SAFE_MODE
#define BLOCK(T) tBlock<T>  // type for a block of a block array - safemode access check
#define PTR_TO_BLOCK(T,M,P,N) tBlock<T>(M,P,N)  // block constructor wrapper - M-block size, P-pointer, N-name
#else
#define BLOCK(X) X*              // type for a block of a block array - release plain pointer
#define PTR_TO_BLOCK(T,M,P,N) P  // block constructor wrapper 
#endif
//----------------------------------------------------------------------------------------------------------------------
// Fixed-size small vector - unallocatable, must be connected to a given pointer, 
// Intended to be used as a return-type in access to block vectors via operator[]
// - fixed size given at constructor. no allocations!
// - explicit access check at safe mode
// - knows name of its creator
//----------------------------------------------------------------------------------------------------------------------
template <typename T> class tBlock {
private:
    T *V; // data array
    const char *name; // pointer to parent object's name.
    int N; // size (32bit only, because it is supposed to be a small vector)

    inline tBlock(){N=0; V=NULL;} // not to be used anywhere  
    inline tBlock<T>& operator=(const tBlock<T>& object){ // forbids assignments from outside
        if(this != &object){ N=object.N; V=object.V; } 
        return *this;
    }

public:
    inline void Reset(){ N=0; V=NULL; name=NULL; }
    inline tBlock(const tBlock<T> &b){ N=b.N; V=b.V; name = b.name; }
    inline tBlock(int n, T *v, const char *lbl=NULL){ N=n; V=v; name = lbl;} 
    inline~tBlock(){}
    
    inline operator       T* () { return V; }
    inline operator const T* () const { return V; }

    inline int Size() const { return N; }
    inline T*  Body() const { return V; }

    inline bool Allocated(void) const { return N>0 && V!=NULL; }

    inline T& operator[](int i){      
        SAFE_ASSERT((V && i>=0 && i<N ), "tBlock: wrong index [%d] %d (%s)", i, N, name ? name : "");  
        return V[i];
    }
    inline const T& operator[](int i)const{ 
        SAFE_ASSERT((V && i>=0 && i<N ), "tBlock: wrong index [%d] %d (%s)", i, N, name ? name : ""); 
        return V[i];
    }
};


//-----------------------------------------------------------------------------------------------------------------------
// Basic 1D vector of data - explicit access check, alloc via mem manager, has a name 
// Intended to be used for big long-living arrays at higher level
// - explicit access check at safe mode
// - fixed size at allocation: new elements cannot be added 
// - alloc via mem manager GimmeMem
// - knows its name for debug reports
//-----------------------------------------------------------------------------------------------------------------------
template <tIdxKind I, typename T> class tArray { // Basic 1D vector of data

private:
    tArray(const tArray<I,T> &b); // force compiler not to create default copy constructor    
    //inline T &operator[](int i){ return VV[i];} // forbids access with plain ints
    //inline const T &operator[](int i) const { return VV[i]; } // forbids access with plain ints

protected:
    T *VV;       // 1D array - vector values.
    size_t NN;   // size 
    string name; // text label (for memory reports and error messages)

    // reset to default all fields
    inline void init_vec(string lbl=""){ NN=0; VV=NULL; name = lbl; }

public:

// Simple getters/setters ---------------------------------------------------------------------------------------------
    inline tInt Size() const { return (tInt)NN; }
    inline T* Body() const { return VV; }
    inline bool Allocated(void)const{ return NN>0; };
    inline void SetName(const char* lbl){ name = lbl ? string(lbl) : ""; }

// Allocations -------------------------------------------------------------------------------------------------------- 
    inline void Alloc(tIdx<I> n, string lbl=""){
        ASSERT(NN==0,"tArray::Alloc %s: already allocated", lbl.c_str()); 
        ASSERT(n>0 , "tArray::Alloc %s: wrong size %llu",   lbl.c_str(), ULL(n.idx())); 
        VV=GimmeMem<T>(n.idx(), lbl.c_str()); NN=(size_t)n.idx(); name = lbl;
    }
    //inline void Alloc(tIdx<I> n, string lbl=""){ Alloc(n.idx(),lbl); }
    inline void Dealloc(){ if(NN>0 && VV) FreeMem(VV); init_vec(""); }

// Basic constructors -------------------------------------------------------------------------------------------------
    tArray(string lbl=""){ init_vec(lbl); }
    tArray(tIdx<I> n, string lbl=""){ init_vec(); Alloc(n, lbl); } 

// Destructor and deallocs --------------------------------------------------------------------------------------------
   ~tArray(){ Dealloc(); }

// Operators ----------------------------------------------------------------------------------------------------------
    inline tArray<I,T>& operator=(tArray<I,T> &object){ // sequential!
        if(this == &object) return *this;
        name = object.name;
        tIdx<I> n; n=NN;
        if(!Allocated() && object.Allocated()) Alloc(n, name.c_str()); //if not yet allocated
        ASSERT(object.NN == NN, "tArray<I,T>: assignment of arrays with different size %llu %llu (%s %s)", 
               ULL(NN), ULL(object.NN), name.c_str(), object.name.c_str());
        for(size_t i=0; i<NN; i++) VV[i]=object.VV[i];
        return *this;
    }
    inline T &operator[](const tIdx<I> &i){        
        SAFE_ASSERT(i>=0 && i<(tInt)NN, "tArray %s[%llu] out of size %llu", name.c_str(), ULL(i.idx()), ULL(NN));
        return VV[i.idx()];
    }
    inline const T &operator[](const tIdx<I> &i) const { // const version
        SAFE_ASSERT(i>=0 && i<(tInt)NN, "tArray %s[%llu] out of size %llu", name.c_str(), ULL(i.idx()), ULL(NN)); 
        return VV[i.idx()];
    } 
};

//----------------------------------------------------------------------------------------------------------------------
// Basic block vector with equal block size for all elements - explicit access check, alloc via mem manager, has name 
// Default data type for block 2D arrays NxM on a base of 1D vector
// Intended to be used for big long-living block arrays at higher level
// - explicit access check at safe mode
// - fixed size at allocation: new elements cannot be added 
// - alloc via mem manager GimmeMem
// - knows its name for debug reports
//----------------------------------------------------------------------------------------------------------------------
template <tIdxKind I, typename T> class tBlockArray : public tArray<I,T>{

private:
    tBlockArray(const tBlockArray<I,T> &b); // force compiler not to create default copy constructor    

protected:
    tInt N; // number of blocks
    int M;  // size of block

    // reset to default all fields
    inline void init_bvec(string lbl=""){ tArray<I,T>::init_vec(lbl); N=0; M=0; } 

    inline void setup_bvec(tInt n, int m, T *v){ // to correctly set sizes
        N=n; M=m; tArray<I,T>::NN=N*M; tArray<I,T>::VV=v; 
    } 

public:

// Simple getters -----------------------------------------------------------------------------------------------------
    inline tInt Size() const {return N;}     // Size of block vector - number of blocks
    inline int BlockSize() const {return M;} // size of each block
    inline size_t BodySize() const {return tArray<I,T>::NN;} // total size of vector (number of scalar values) 
    inline       tArray<I,T>& BodyVec()       {return ((tArray<I,T>&)(*this)); } // This is to access block vector as vector 
    inline const tArray<I,T>& BodyVec() const {return ((tArray<I,T>&)(*this)); } // to avoid nested loops (do not use Body() for that!) 
    inline bool Allocated()const{return (M>0 && N>0); }

// Allocations -------------------------------------------------------------------------------------------------------- 
    inline void Alloc(tIdx<I> n, int m, string lbl=""){
        ASSERT((n>0 && m>0), "tBlockArray Alloc: wrong input %d %d (%s)", n, m, lbl.c_str());
        ASSERT(!Allocated(), "tBlockArray Alloc: already allocated! (%s)", lbl.c_str());
        size_t Size = (size_t)n.idx() * (size_t)m; 
        setup_bvec(n.idx(), m, GimmeMem<T>(Size, lbl.c_str()));
        tArray<I,T>::name = lbl;
    }
    inline void Dealloc(){ tArray<I,T>::Dealloc(); init_bvec(); }

// Basic constructors -------------------------------------------------------------------------------------------------
    tBlockArray(string lbl=""){ init_bvec(lbl); }
//  tBlockArray(tInt n, int m, string lbl=""){ init_bvec(); Alloc(n, m, lbl); }
    tBlockArray(tIdx<I> n, int m, string lbl=""){ init_bvec(); Alloc(n, m, lbl); }
        
// Destructor and deallocs --------------------------------------------------------------------------------------------
    ~tBlockArray(){ Dealloc(); }
        
// Operators ----------------------------------------------------------------------------------------------------------
    inline tBlockArray<I,T>& operator=(const tBlockArray<I,T>& object){ // sequential
        if(this == &object) return *this;
        tArray<I,T>::name = object.tArray<I,T>::name;
        if(!Allocated() && object.Allocated()) Alloc(object.N, object.M, tArray<I,T>::name.c_str()); //if not yet allocated
        ASSERT(object.N==N, "tBlockArray: assignment of arrays with different size (%s)", 
                            object.N, N, tArray<I,T>::name.c_str());
        ASSERT(object.M==M, "tBlockArray: assignment of arrays with different size of block %d %d (%s)", 
                            object.M, M, tArray<I,T>::name.c_str());
        for(size_t i=0; i<(size_t)N*M; i++) tArray<I,T>::VV[i]=object.VV[i];
        return *this;
    }

    inline BLOCK(T) operator[](tIdx<I> i){      
        SAFE_ASSERT((i>=0 && i<N), "tBlockArray %s[%llu] out of size %llu", tArray<I,T>::name.c_str(), ULL(i.idx()), ULL(N));
        SAFE_ASSERT((N>0&&M>0&&tArray<I,T>::VV), "tBlockArray access to empty array! (%s)", tArray<I,T>::name.c_str());
        return PTR_TO_BLOCK(T, M, (tArray<I,T>::VV+(size_t)(i.idx()*M)), (tArray<I,T>::name.c_str()));    
    }
    inline const tBlock<T> operator[](tIdx<I> i)const {// версия для const & аргумента
        SAFE_ASSERT((i>=0&&i<N), "tBlockArray %s[%llu] out of size %llu", tArray<I,T>::name.c_str(), ULL(i.idx()), ULL(N));
        SAFE_ASSERT((N>0&&M>0&&tArray<I,T>::VV), "tBlockArray access to empty array! (%s)", tArray<I,T>::name.c_str());
        return PTR_TO_BLOCK(T, M, (tArray<I,T>::VV+(size_t)(i.idx()*M)), (tArray<I,T>::name.c_str()));    
    } 

};


//----------------------------------------------------------------------------------------------------------------------
// Basic block vector with variable block size (fixed at construction, not dynamic)
// Intended to be used for CSR topology IA/JA at higher level
// - explicit access check at safe mode
// - fixed size at allocation: new elements cannot be added 
// - alloc via mem manager GimmeMem
// - knows its name for debug reports
//----------------------------------------------------------------------------------------------------------------------
template <tIdxKind I, typename T> class tVBlockArray : public tArray<I,T>{
private:
    tVBlockArray(const tVBlockArray<I,T> &b); // force compiler not to create default copy constructor  

protected:
    tInt N; // size - number of blocks (rows)
    tInt *IA; // row pointer addresses for the beginning of the element data
    // IA - row pointer, tArray::NN = IA[N] 
    // JA - tArray::VV - column pointer
    int sorted; // flag that blocks in V are sorted ascending - to use fast search of position in block 

    inline void init_dbvec(string lbl=""){ tArray<I,T>::init_vec(lbl); N = 0; IA = NULL; sorted = 0; }

    inline void setup_dbvec(tInt n, tInt *ia, T *ja){ // to correctly set sizes
        N=n; IA=ia; tArray<I,T>::NN=ia[N]; tArray<I,T>::VV=ja; 
    } 

public:

// Getters/setters ---------------------------------------------------------------------------------------------
    inline bool Allocated()const { return (N>0); }
    inline tInt Size() const {return N;} // size - number of blocks (rows)
    inline tInt BodySize() const { return (IA ? IA[N] : 0); } // total size of vector in scalar values 
    inline bool IsSorted() const { return sorted>0; } 
    inline void SetSorted(){   sorted=1; }
    inline void SetUnsorted(){ sorted=0; }
    inline       tArray<I,T>& BodyVec()       {return ((tArray<I,T>&)(*this)); } // This is to access block vector as vector 
    inline const tArray<I,T>& BodyVec() const {return ((tArray<I,T>&)(*this)); } // to avoid nested loops (do not use Body() for that!) 
    
    inline int BlockSize(const tIdx<I> i)const{ // blocksize. number of values in i-th block. to be renamed!
        SAFE_ASSERT((i>=0&&i<N), "tVBlockArray GetM: %d out of size %d! (%s)", i.idx(), N, tArray<I,T>::name.c_str());   
        SAFE_ASSERT((N>0&&IA&&tArray<I,T>::VV), "tVBlockArray GetM: access to empty array! (%s)", tArray<I,T>::name.c_str());
        SAFE_ASSERT((IA[i.idx()+1]-IA[i.idx()]>=0 && IA[i.idx()]>=0), "tVBlockArray GetM: wrong IA! (%s)", tArray<I,T>::name.c_str());
        return (int)(IA[i.idx()+1]-IA[i.idx()]); 
    }

    // pointers to data - CSR synonims
    inline       tInt* GetIA()       { SAFE_ASSERT((IA!=NULL), "tVBlockArray GetIA: IA is NULL!(%s)", tArray<I,T>::name.c_str()); return IA;}
    inline const tInt* GetIA() const { SAFE_ASSERT((IA!=NULL), "tVBlockArray GetIA: IA is NULL!(%s)", tArray<I,T>::name.c_str()); return IA;}
    inline       T* GetJA()       { return tArray<I,T>::Body();} 
    inline const T* GetJA() const { return tArray<I,T>::Body();} 

    // finds position of j-th column in i-th row
    inline int FindIndex(const tIdx<I> i, T j, int doCrash=1) const { // returns position in 1D vector (returns JA)
        SAFE_ASSERT((i>=0&&i<N), "tVBlockArray FindIndex: i=%d out of size %d! (%s)", i, N, tArray<I,T>::name.c_str());
        int res =  sorted ? findCSRIndex(IA, tArray<I,T>::VV, i, j) :
                            findCSRIndexSucc(IA, tArray<I,T>::VV, i, j); 
        ASSERT((res>=0 || !doCrash), "tVBlockArray FindIndex: in row %d column %d not found (%s)", i, j, tArray<I,T>::name.c_str());
        return res; 
    }    
    inline int FindIndexPos(tIdx<I> i, T j, int doCrash=1) const { // returns position in i-th block
        return FindIndex(i, j, doCrash)-IA[i]; 
    }

 // Basic constructors -------------------------------------------------------------------------------------------------
    tVBlockArray(string lbl=""){ init_dbvec(lbl); }
    tVBlockArray(tInt n, int m, string lbl=""){ init_dbvec(); Alloc(n, m, lbl); }
   ~tVBlockArray(){ Dealloc(); } 

   inline void Dealloc(){ 
        if(tArray<I,T>::VV) FreeMem(tArray<I,T>::VV); 
        if(IA)FreeMem(IA); 
        init_dbvec("");
    }

// Allocations -------------------------------------------------------------------------------------------------------- 
    inline void Alloc(tIdx<I> n, int m, string lbl="") {
        if(Allocated())   crash("tVBlockArray Alloc: already allocated! (%s)", lbl.c_str());
        if((n<=0)||(m<=0))crash("tVBlockArray Alloc: wrong input %d %d (%s)", n, m, lbl.c_str());
        size_t nd = (size_t)n.idx() * (size_t)m; 
        if(sizeof(tInt)==4) if(nd>=MAX_SIGNED_INT) // check that body size doesn't exceed int range - because IA is int
            crash("tVBlockArray Alloc: body size exceeds 32bit range (%s)", lbl.c_str());        

        N = n.idx();
        IA = GimmeMem<int>(N+1, lbl.c_str());
        for(int i=0; i<N+1; i++) IA[i] = i*m;

        setup_dbvec(n.idx(), IA, GimmeMem<T>(nd, lbl.c_str()) ); 
        sorted = 0; 
        tArray<I,T>::name = lbl;
    }

    template <typename idxt> inline void Alloc(tIdx<I> n, const tArray<I,idxt> &blockSizes, string lbl="") { //idxt - can be char to save space
        ASSERT(!Allocated(), "tVBlockArray Alloc: already allocated! (%s)", lbl.c_str());
        ASSERT(n>0 && blockSizes.Size()==n, "tVBlockArray Alloc: wrong input %d (%s)", n, lbl.c_str());
        size_t nd = 0;
        IA = GimmeMem<tInt>(n.idx()+1, lbl.c_str());
        IA[0] = 0;
        tIdx<I> i;
        for(i=0; i<n.idx(); i++){ 
            if((tInt)blockSizes[i]<0)crash("tVBlockArray Alloc: wrong block size %d (%s)", blockSizes[i], lbl.c_str());
            nd += blockSizes[i];
            IA[i.idx()+1] = (tInt)nd;
        }
        if(nd==0) crash("tVBlockArray Alloc: empty body (%s)", lbl.c_str());
        if(sizeof(tInt)==4) if(nd>=MAX_SIGNED_INT) // check that body size doesn't exceed int range - because IA is int
            crash("tVBlockArray Alloc: body size exceeds 32bit range (%s)", lbl.c_str());    

        setup_dbvec(n.idx(), IA, GimmeMem<T>(nd, lbl.c_str()) );
        sorted = 0; 
        tArray<I,T>::name = lbl;
    }
 
//// Construction based on external data --------------------------------------------------------------------------------
//// Data pointer v will belong to this object, so it must not be deallocated elsewhere!
//// Data pointer v WILL BE deallocated at Dealloc or destructor (so it mustn't be used elsewhere) 
//    inline void Construct(tInt n, tInt *ia, T *v/*can be null*/, string lbl="") {
//        if(v==NULL && n>0 && ia!=NULL) if(ia[n]>0) v = GimmeMem<T>(ia[n], lbl.c_str());
//
//        if(Allocated())crash("tVBlockArray Construct: already allocated! (%s)", lbl.c_str());
//        if(n<=0||!ia)  crash("tVBlockArray Construct: wrong input %d (%s)", n, lbl.c_str());
//        if(ia[0]!=0)   crash("tVBlockArray Construct: wrong IA[0]!=0 (%s)", lbl.c_str());
////      if(ia[n]<=0)   crash("tVBlockArray Alias: wrong IA[N]<=0 (%s)", lbl.c_str());
//        if(v==NULL && ia[n]>0) crash("tVBlockArray Construct: v is null! (%s)", n, lbl.c_str());
//        for(tInt i=0; i<n; i++) if(ia[i+1]-ia[i]<0 || ia[i+1]<0/*to avoid runout of int range*/) 
//            crash("tVBlockArray Construct: wrong IA (%s)", lbl.c_str());
//
//        setup_dbvec(n, ia, v);
//        sorted = 0;
//        tArray<I,T>::name = lbl;
//    }

// Operators ----------------------------------------------------------------------------------------------------------

    inline tVBlockArray<I,T>& operator=(const tVBlockArray<I,T>& object){
        if(this == &object) return *this;
        if(Allocated()) Dealloc(); 
        if(object.Allocated()) ConstructCopy(object.N, object.IA, object.VV, object.name.c_str()); 
        sorted = object.sorted;
        return *this;
    }
    inline tVBlockArray<I,T>& operator=(T x){
        if(!Allocated()) crash("tVBlockArray =: not allocated! (%s)", tArray<I,T>::name.c_str());
        for(int i=0; i<BodySize(); i++) tArray<I,T>::VV[i]=x;
        return *this;
    }

    inline BLOCK(T) operator[](tIdx<I> i){     
        SAFE_ASSERT((i<MAX_SIGNED_INT), "DtBlockArray [] out of int range! (%s)", tArray<I,T>::name.c_str());
        SAFE_ASSERT((i>=0&&i<N), "tVBlockArray [%d] out of size %d! (%s)", i, N, tArray<I,T>::name.c_str()); 
        SAFE_ASSERT((N>0&&IA&&tArray<I,T>::VV), "tVBlockArray access to empty array! (%s)", tArray<I,T>::name.c_str());
        return PTR_TO_BLOCK(T, BlockSize(i), (tArray<I,T>::VV+IA[i.idx()]), (tArray<I,T>::name.c_str()));  
    }
    inline const BLOCK(T) operator[](tIdx<I> i)const { // версия для const & аргумента
        SAFE_ASSERT((i<MAX_SIGNED_INT), "tVBlockArray [] out of int range! (%s)", tArray<I,T>::name.c_str());
        SAFE_ASSERT((i>=0&&i<=N), "tVBlockArray [%d] out of size %d! (%s)", i, N, tArray<I,T>::name.c_str()); 
        SAFE_ASSERT((N>0&&IA&&tArray<I,T>::VV), "tVBlockArray access to empty array! (%s)", tArray<I,T>::name.c_str());
        return PTR_TO_BLOCK(T, BlockSize(i), (tArray<I,T>::VV+IA[i.idx()]), (tArray<I,T>::name.c_str()));  
    }
};


void DoSomething1(tIdx<IDX_N> /*in*/, BLOCK(double) /*v*/){} 
void DoSomething2(tIdx<IDX_N> /*in*/, tIdx<IDX_E> /*ie*/, BLOCK(double) /*v*/){}


tInt main(tInt, const char**){
#if(0)
    tIdx<IDX_N> in; in = 0;
    tIdx<IDX_E> ie; ie = 0;

    in = in + 10;
//  in = 10 + in; // error
    ie = (ie + 1)*2;
    in = in + in; 
//  in = in + ie; // error
    
    tIdx<IDX_N> n0, nn; n0=0; nn = 100;
    tIdx<IDX_E> e0, en; e0=0; en = 100;
    tArray<IDX_N, int> AN(nn, "AN");
    tArray<IDX_E, int> AE(en, "AE");
        
    for(tIdx<IDX_N> i=n0; i<nn; i++) AN[i]=0;
    for(tIdx<IDX_E> i=e0; i<en; i++) AE[i]=0;
//  for(tIdx<IDX_E> i=0; i<100; i++) AN[i]=0; // error
//  for(tInt        i=0; i<100; i++) AN[i]=0; // error
//  for(tIdx<IDX_N> i=0; i<101; i++) AN[i]=0; // error (runtime)
    
//  AN[ 0] = 0; // error
//  AN[ie] = 0; // error
    AN[in] = 0;
//  AE[in] = 0; // error
    AE[ie] = 0;

    tArray<IDX_N, int> AN2;
    AN2 = AN;
//  AN2 = AE; // error 
    

    tBlockArray<IDX_N, tIdx<IDX_E> > BANE; 
    BANE.Alloc(nn, 5, "BANE");
    tBlockArray<IDX_E, tIdx<IDX_N> > BAEN(en, 5, "BAEN");

    for(tIdx<IDX_N> i=n0; i<nn; i++){
        for(int j=0; j<5; j++) BANE[i][j] = 0;
//      for(int j=0; j<5; j++) BANE[i][j] = in; // error
        for(int j=0; j<6; j++) BANE[i][j] = 0;  // error (runtime)
        for(int j=0; j<5; j++) BANE[i][j] = ie; 
    }

//  for(tIdx<IDX_E> i=0; i<100; i++) // error
//      for(int j=0; j<5; j++) BANE[i][j] = 0;

    AN[BAEN[ie][0]]=0;
//  AN[BANE[ie][0]]=0; // error
//  AN[BANE[in][0]]=0; // error
//  AN[BANE[in][7]]=0; // error (runtime); 
        

    tVBlockArray<IDX_E, tIdx<IDX_N> > VBAEN(100, 4, "VBAEN");
    tVBlockArray<IDX_N, tIdx<IDX_E> > VBANE;

    tArray<IDX_N,char> bvsizes(nn,"bvsizes"); 
    for(tIdx<IDX_N> i=n0; i<nn; i++) bvsizes[i] = (char)(i.idx()%4); 
    VBANE.Alloc(100,bvsizes,"VBANE"); 
//  VBANE.Alloc(101,bvsizes,"VBANE"); // error (runtime)
//  VBAEN.Alloc(100,bvsizes,"VBAEN"); // error

//  for(tIdx<IDX_N> i=0; i<100; i++) // error
//      for(int j=0; j<4; j++) VBAEN[i][j] = 0;
    
    for(tIdx<IDX_E> i=e0; i<en; i++){
        for(int j=0; j<4; j++) VBAEN[i][j] = 0;
//      for(int j=0; j<4; j++) VBAEN[i][j] = ie; // error
    }
    for(tIdx<IDX_N> i=n0; i<nn; i++){
        for(int j=0; j<VBANE.BlockSize(i); j++) VBANE[i][j] = 0;
//      for(int j=0; j<4; j++) VBAEN[i][j] = ie; // error
    }
    

    ie = 3; 
    in = 2; 
    AN[VBAEN[ie][0]]=0;
//  AN[VBANE[ie][0]]=0; // error
//  AN[VBANE[in][0]]=0; // error
//  AN[VBANE[in][7]]=0; // error (runtime); 
    
    BANE[VBAEN[ie][0]][0]=0;
//  BANE[VBANE[in][0]][0]=0; // error

    VBANE[BAEN[ie][0]][0]=0;
//  VBANE[BANE[in][0]]=0; // error

    printf("index N %d, index E %d\n", (int)in.idx(), (int)ie.idx()); 
    printf("index %d\n", (int)BANE.Size()); 
    printf("index %d\n", (int)BAEN.Size()); 
    printf("index %d\n", (int)VBANE.Size()); 
    printf("index %d\n", (int)VBAEN.Size()); 


    tIdx<IDX_N> nN=nn; // number of mesh nodes
    tIdx<IDX_E> nE=en; // number of mesh elements    
    // ...
    // adjacency/topology data
    tVBlockArray<IDX_N, tIdx<IDX_E> > NE_topo("NE_topo"); // for each node stores elements it belings to
    tVBlockArray<IDX_E, tIdx<IDX_N> > EN_topo("EN_topo"); // for each element stores its nodes
    // ...
    tArray<IDX_E,double> Evol(nE, "Evol"); // volumes of mesh elements
    tBlockArray<IDX_N,double> Ncoor(nN, 3, "Nvol"); // coordinates of mesh nodes
    //tBlockArray<IDX_N,double> Ncoor(nE, "Nvol"); // error
    // ...
    void DoSomething(tIdx<IDX_E> i, BLOCK(double) c); 

 // for(tIdx<IDX_N> i=0; i<nN; i++){ // error
    for(tIdx<IDX_E> i=e0; i<nE; i++){
//      for(int j=0; j<NE_topo.BlockSize(i); j++){ //error
        for(int j=0; j<EN_topo.BlockSize(i); j++){
//          DoSomething(i,Ncoor[NE_topo[i][j]]); //error
//          DoSomething(i,Ncoor[EN_topo[j][i]]); //error 
            DoSomething(i,Ncoor[EN_topo[i][j]]); 
        }
    }

#else 

  tIdx<IDX_N> nn; // number of nodes in MPI extended subdomain 
  tIdx<IDX_N> n_beg, n_end; // range for extended subdomain (owned+halo nodes) 
  nn = 1000; 
  n_beg=0; n_end = nn; 
  // subsets (vertex-centred case, mesh decomposition for the mesh nodal graph)
  //tIdx<IDX_N> n_owned_beg, n_owned_end; // range for owned nodes (inner+interface)
  //tIdx<IDX_N> n_inner_beg, n_inner_end; // range for inner nodes
  //tIdx<IDX_N> n_iface_beg, n_iface_end; // range for interface nodes
  //tIdx<IDX_N> n_halo_beg,  n_halo_end;  // range for halo nodes
  // ...
  tIdx<IDX_E> en; // number of elements in MPI subdomain 
  tIdx<IDX_E> e_beg, e_end; // range for elements
  en = 1000; 
  e_beg=0; e_end = en; 
                            // ...
  tBlockArray<IDX_N, double> BAN("BAN"); // some block array over nodes 
  tVBlockArray<IDX_E, tIdx<IDX_N> > EN_topo; // EN topo - a variable-block array over elements

  BAN.Alloc(nn,3,"BAN"); 
  EN_topo.Alloc(en,4,"EN_topo"); 
  for(tIdx<IDX_E> ie=e_beg; ie<e_end; ++ie){  // loop over elements
      for(int j=0; j<EN_topo.BlockSize(ie); ++j){ 
          EN_topo[ie][j] = (ie.idx()+j) % nn.idx();
      }
  }
  // ...

  // ...
//for(tIdx<IDX_N> in=e_beg; in<n_end; ++ie) // compile-time error
//for(tIdx<IDX_E> in=n_beg; in<n_end; ++ie) // compile-time error
  for(tIdx<IDX_N> in=n_beg; in<n_end; ++in) // loop over nodes
      DoSomething1(in, BAN[in]);
  // ...
  for(tIdx<IDX_E> ie=e_beg; ie<e_end; ++ie){  // loop over elements
      tIdx<IDX_N> in;
      BLOCK(tIdx<IDX_N>) nodes = EN_topo[ie]; // the nodes of the element
//    BLOCK(tIdx<IDX_N>) nodes = EN_topo[in]; // compile-time error
//    BLOCK(tIdx<IDX_E>) nodes = EN_topo[ie]; // compile-time error
      for(int j=0; j<EN_topo.BlockSize(ie); ++j){ 
//    for(int j=0; j<=EN_topo.BlockSize(ie); ++j){ // runtime error at block access check 
          in = nodes[j];
//        const tIdx<IDX_E> in = nodes[j]; // compile-time error
          DoSomething2(in, ie, BAN[in]); 
//        DoSomething2(ie, in, BAN[in]); // compile-time error
      }
    }  
     
    tIdx<IDX_N> NNN_BEG, NNN_END; NNN_BEG=0; NNN_END=nn; 
    //tIdx<IDX_N> myNNN = NNN_BEG; 
    const int MMM = 8; 
    double buf = 0.0;
    tBlockArray<IDX_N, double > TEST1(NNN_END, MMM, "TEST1"); 
    tBlockArray<IDX_N, double > TEST2(NNN_END, MMM, "TEST1"); 

    tArray<IDX_C, double> A; // Array over cells
    tIdx<IDX_N> in; 
    tIdx<IDX_E> ie; 
    tIdx<IDX_S> is; 
    
    A.Alloc(nn,"A");
    in=n_beg;
    A[in] = 0.0; 
//  A[ie] = 0.0; // Runtime error
//  A[is] = 0.0; // Link-time error
    //tBlockArray<IDX_S, double > TEST4(NNN_END, MMM, "TEST1"); 

    //CR_Timer.cr_reset();
    //CR_Timer.cr_beg("TEST1");
    for(int k=0; k<100; k++){
        for(tIdx<IDX_N> i=NNN_BEG; i<NNN_END; ++i){
            BLOCK(double) B1 = TEST1[i];
            BLOCK(double) B2 = TEST2[i];
            for(int j=0; j<MMM; ++j){ B1[j] = 3*i.idx(); B2[j] = 2*i.idx(); }
        }
        for(tIdx<IDX_N> i=NNN_BEG; i<NNN_END; ++i){
            BLOCK(double) B1 = TEST1[i];
            BLOCK(double) B2 = TEST2[i];
            for(int j=0; j<MMM; ++j) B1[j] += B2[j]; 
        }
        for(tIdx<IDX_N> i=NNN_BEG; i<NNN_END; ++i){
            BLOCK(double) B1 = TEST1[i];
            for(int j=0; j<MMM; ++j) buf += B1[j]; 
        }
    }
    //CR_Timer.cr_end("TEST1");

    TEST1.Dealloc();
    TEST2.Dealloc();

    const int NNN = nn.idx();
    double *DTEST1 = new double[NNN*MMM]; 
    double *DTEST2 = new double[NNN*MMM]; 

    //CR_Timer.cr_beg("TEST2");
    for(int k=0; k<100; k++){
        for(int i=0; i<NNN; ++i){
            double *B1 = DTEST1 + i*MMM;
            double *B2 = DTEST2 + i*MMM;
            for(int j=0; j<MMM; ++j){ B1[j] = 3*i; B2[j] = 2*i; }
        }
        for(int i=0; i<NNN; ++i){
            double *B1 = DTEST1 + i*MMM;
            double *B2 = DTEST2 + i*MMM;
            for(int j=0; j<MMM; ++j) B1[j] += B2[j]; 
        }
        for(int i=0; i<NNN; ++i){
            double *B1 = DTEST1 + i*MMM;
            for(int j=0; j<MMM; ++j) buf += B1[j]; 
        }
    }
    //CR_Timer.cr_end("TEST2");

    TEST1.Dealloc();

    //CR_Timer.cr_info_short(1);
    printf("bufff: %g\n", buf);
#endif

}
