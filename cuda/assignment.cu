#define SQR(A) ((A)*(A))

__device__ float caf(int i, float x, int cao_value) {
  switch (cao_value) {
  case 1 : return 1.0;
  case 2 : {
    switch (i) {
    case 0: return 0.5-x;
    case 1: return 0.5+x;
    default:
      return 0.0;
    }
  } 
  case 3 : { 
    switch (i) {
    case 0: return 0.5*SQR(0.5 - x);
    case 1: return 0.75 - SQR(x);
    case 2: return 0.5*SQR(0.5 + x);
    default:
      return 0.0;
    }
  case 4 : { 
    switch (i) {
    case 0: return ( 1.0+x*( -6.0+x*( 12.0-x* 8.0)))/48.0;
    case 1: return (23.0+x*(-30.0+x*(-12.0+x*24.0)))/48.0;
    case 2: return (23.0+x*( 30.0+x*(-12.0-x*24.0)))/48.0;
    case 3: return ( 1.0+x*(  6.0+x*( 12.0+x* 8.0)))/48.0;
    default:
      return 0.0;
    }
  }
  case 5 : {
    switch (i) {
    case 0: return (  1.0+x*( -8.0+x*(  24.0+x*(-32.0+x*16.0))))/384.0;
    case 1: return ( 19.0+x*(-44.0+x*(  24.0+x*( 16.0-x*16.0))))/ 96.0;
    case 2: return (115.0+x*       x*(-120.0+x*       x*48.0))  /192.0;
    case 3: return ( 19.0+x*( 44.0+x*(  24.0+x*(-16.0-x*16.0))))/ 96.0;
    case 4: return (  1.0+x*(  8.0+x*(  24.0+x*( 32.0+x*16.0))))/384.0;
    default:
      return 0.0;
    }
  }
  case 6 : {
    switch (i) {
    case 0: return (  1.0+x*( -10.0+x*(  40.0+x*( -80.0+x*(  80.0-x* 32.0)))))/3840.0;
    case 1: return (237.0+x*(-750.0+x*( 840.0+x*(-240.0+x*(-240.0+x*160.0)))))/3840.0;
    case 2: return (841.0+x*(-770.0+x*(-440.0+x*( 560.0+x*(  80.0-x*160.0)))))/1920.0;
    case 3: return (841.0+x*(+770.0+x*(-440.0+x*(-560.0+x*(  80.0+x*160.0)))))/1920.0;
    case 4: return (237.0+x*( 750.0+x*( 840.0+x*( 240.0+x*(-240.0-x*160.0)))))/3840.0;
    case 5: return (  1.0+x*(  10.0+x*(  40.0+x*(  80.0+x*(  80.0+x* 32.0)))))/3840.0;
    default:
      return 0.0;
    }
  }
  case 7 : {
    switch (i) {
    case 0: return (    1.0+x*(   -12.0+x*(   60.0+x*( -160.0+x*(  240.0+x*(-192.0+x* 64.0))))))/46080.0;
    case 1: return (  361.0+x*( -1416.0+x*( 2220.0+x*(-1600.0+x*(  240.0+x*( 384.0-x*192.0))))))/23040.0;
    case 2: return (10543.0+x*(-17340.0+x*( 4740.0+x*( 6880.0+x*(-4080.0+x*(-960.0+x*960.0))))))/46080.0;
    case 3: return ( 5887.0+x*          x*(-4620.0+x*         x*( 1680.0-x*        x*320.0)))   /11520.0;
    case 4: return (10543.0+x*( 17340.0+x*( 4740.0+x*(-6880.0+x*(-4080.0+x*( 960.0+x*960.0))))))/46080.0;
    case 5: return (  361.0+x*(  1416.0+x*( 2220.0+x*( 1600.0+x*(  240.0+x*(-384.0-x*192.0))))))/23040.0;
    case 6: return (    1.0+x*(    12.0+x*(   60.0+x*(  160.0+x*(  240.0+x*( 192.0+x* 64.0))))))/46080.0;
    default:
      return 0.0;
    }
  }
  }}
  return 0.0;
}

__device__ int wrap_index(int ind, const int mesh) {
  if(ind < 0)
    return ind + mesh;
  else if(ind >= mesh)
    return ind - mesh;
  else 
    return ind;	   
}

__device__ inline void atomicadd(float* address, float value){
#if !defined __CUDA_ARCH__ || __CUDA_ARCH__ >= 200 // for Fermi, atomicAdd supports floats
  atomicAdd(address, value);
#elif __CUDA_ARCH__ >= 110
#warning Using slower atomicAdd emulation
  // float-atomic-add from 
  // [url="http://forums.nvidia.com/index.php?showtopic=158039&view=findpost&p=991561"]
  float old = value;
  while ((old = atomicExch(address, atomicExch(address, 0.0f)+old))!=0.0f);
#else
#warning I need at least compute capability 1.1
#endif
}

#define r_ind(A,B,C) ((A)*m_size*m_size + (B)*m_size + (C))

__global__ void assign_charges(const float * const pos, const float * const q,
float *mesh, const int m_size, const int cao, const float pos_shift, const
float hi) {
      /** id of the particle **/
      int id = blockIdx.x;
      /** position relative to the closest gird point **/
      float m_pos[3];
      /** index of the nearest mesh point **/
      int nmp_x, nmp_y, nmp_z;      

      m_pos[0] = pos[3*id + 0] * hi - pos_shift;
      m_pos[1] = pos[3*id + 1] * hi - pos_shift;
      m_pos[2] = pos[3*id + 2] * hi - pos_shift;

      nmp_x = (int) floor(m_pos[0] + 0.5);
      nmp_y = (int) floor(m_pos[1] + 0.5);
      nmp_z = (int) floor(m_pos[2] + 0.5);

      m_pos[0] -= nmp_x;
      m_pos[1] -= nmp_y;
      m_pos[2] -= nmp_z;

      nmp_x = wrap_index(nmp_x + threadIdx.x, m_size);
      nmp_y = wrap_index(nmp_y + threadIdx.y, m_size);
      nmp_z = wrap_index(nmp_z + threadIdx.z, m_size);
      /** should be atomic **/
      mesh[r_ind(nmp_x, nmp_y, nmp_z)].x += caf(threadIdx.x, m_pos[0], cao)*caf(threadIdx.y, m_pos[1], cao)*caf(threadIdx.z, m_pos[2], cao)*q[id];
}

