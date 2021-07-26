#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <math.h>
#include <windows.h>
#include <iostream>
using namespace std;

_Dcomplex operator+(_Dcomplex a, _Dcomplex b) {
 _Dcomplex tmp;
 tmp._Val[0] = a._Val[0] + b._Val[0];
 tmp._Val[1] = a._Val[1] + b._Val[1];
 return tmp;
}
_Dcomplex operator-(_Dcomplex a, _Dcomplex b) {
 _Dcomplex tmp;
 tmp._Val[0] = a._Val[0] - b._Val[0];
 tmp._Val[1] = a._Val[1] - b._Val[1];
 return tmp;
}
_Dcomplex operator*(_Dcomplex a, _Dcomplex b) {
 _Dcomplex tmp;
 tmp._Val[0] = a._Val[0] * b._Val[0] - a._Val[1] * b._Val[1];
 tmp._Val[1] = a._Val[1] * b._Val[0] + a._Val[0] * b._Val[1];
 return tmp;
}
_Dcomplex operator*(_Dcomplex a, double b) {
 _Dcomplex tmp;
 tmp._Val[0] = a._Val[0] * b;
 tmp._Val[1] = a._Val[1] * b;
 return tmp;
}
class fftmath {
protected:
 const _Dcomplex mpii = { 0,-3.14159265358979323846 };
 int outsize = 0, n = 0;
 int overlap;
 double *lfreq, *rfreq, *centerfreq;
 _Dcomplex *out = 0, *cexpb = 0;
 double* dout = 0;
 fftmath(int overlap, int n) :overlap(overlap), n(n) {
  if (cexpb) free(cexpb);
  cexpb = (_Dcomplex*)malloc(sizeof(_Dcomplex)*n);
  for (int i = 0; i < n; i++) cexpb[i] = cexp(mpii * (static_cast<double>(i) / static_cast<double>(n)));
 }
 int log2n(int n) {
  int ret = -1;
  while (n) ret++, n >>= 1;
  return ret;
 }
 virtual _Dcomplex* ft(_Dcomplex *in) = 0;
 virtual double* operator() (_Dcomplex*in) = 0;
 virtual ~fftmath() {
  if (cexpb) free(cexpb);
 }
public:
 double* getlfreq() {
  return lfreq;
 }
 double* getrfreq() {
  return rfreq;
 }
 double* getcenterfreq() {
  return centerfreq;
 }
 int getoutsize() {
  return outsize;
 }
 int getinsize() {
  return n;
 }
};
class lftcalc :virtual public fftmath {
private:
 int *ovldata = 0;
 int **lyf;
 int **lyr;
 int *nly;
 int nnly;
 int _step;
 int nnly2;
 _Dcomplex *tmp1, *tmp2, *tmp3, *sw;
 _Dcomplex* mem1, *mem2;
 double **que;
 int qptr;
 int  r, s;
 double *dovldata;
public:
 lftcalc(int overlap, int n, int r, int s) : r(r), s(s), fftmath(overlap, n) {
  if (overlap < 1) throw exception("overlap must be positive integer");
  if (n&(n - 1)) throw exception("n must be power of 2");
  if (n < 2) throw exception("n must be larger than 1");
  if (r&(r - 1)) throw exception("r must be power of 2");
  if (r < 2) throw exception("r must be larger than 1");
  if (s&(s - 1)) throw exception("s must be power of 2");
  if (s < 1) throw exception("s must be larger than 0");
  if (s >= r) throw exception("s must be smaller than r/2");
  nnly = log2n(n);
  lyf = (int**)malloc(sizeof(int*)*nnly);
  lyr = (int**)malloc(sizeof(int*)*nnly);
  nly = (int*)malloc(sizeof(int)*nnly);
  for (int i = 0; i < nnly; i++) {
   nly[i] = nnly - i;
   lyf[i] = (int*)malloc(sizeof(int)*nly[i]);
   lyr[i] = (int*)malloc(sizeof(int)*nly[i]);
  }
  mem1 = (_Dcomplex*)malloc(sizeof(_Dcomplex)*n);
  mem2 = (_Dcomplex*)malloc(sizeof(_Dcomplex)*n);
  int t = 0;
  int st = r / s;
  _step = log2n(st);
  lyf[0][nly[0] - 1] = 0;
  lyr[0][nly[0] - 1] = r;
  outsize = r;
  int q;
  t += _step;
  while (1) {
   lyf[t][nly[t] - 1] = s;
   q = r << t;
   if (q >= (n >> 1)) {
    lyr[t][nly[t] - 1] = r * (n >> 1) / q;
    outsize += lyr[t][nly[t] - 1] - s;
    nnly2 = t;
    break;
   }
   else lyr[t][nly[t] - 1] = r;
   t += _step;
   outsize += r - s;
  }
  for (int i = 0; i < nnly; i++) {
   for (int j = nly[i] - 2; j >= 0; j--) {
    if (lyr[i][j + 1] > (n >> (i + 1))) {
     lyf[i][j] = 0;
     lyr[i][j] = n >> i;
    }
    else {
     lyf[i][j] = lyf[i][j + 1] << 1;
     lyr[i][j] = lyr[i][j + 1] << 1;
    }
   }
  }
  out = (_Dcomplex*)malloc(sizeof(_Dcomplex)*outsize);
  que = (double**)malloc(sizeof(double*)*overlap);
  for (int i = 0; i < overlap; i++) que[i] = (double*)malloc(sizeof(double)*outsize);
  dout = (double*)malloc(sizeof(double)*outsize);
  lfreq = (double*)malloc(sizeof(double)*outsize);
  rfreq = (double*)malloc(sizeof(double)*outsize);
  centerfreq = (double*)malloc(sizeof(double)*outsize);
  ovldata = (int*)malloc(sizeof(int)*outsize);
  dovldata = (double*)malloc(sizeof(double)*outsize);
  t = 0; int u = 0;
  int v = r / s;
  int w = 0, y = 1, z;
  while (1) {
   for (int a = w; a < r; a++, u += y, t++) {
    if (t == outsize) goto T;
    lfreq[t] = (double)u / (double)(n >> 1);
    rfreq[t] = lfreq[t] + y / (double)(n >> 1);
    centerfreq[t] = (lfreq[t] + rfreq[t])*0.5;
    z = y > overlap ? overlap : y;
    ovldata[t] = overlap / z;
    dovldata[t] = sqrt((double)y)/(double)z;
   }
   y *= v;
   w = s;
  }
 T:
  for (int i = 0; i < overlap; i++)
   for (int j = 0; j < outsize; j++) que[i][j] = 0.;
 }
 _Dcomplex* ft(_Dcomplex *in) {
  int step, tstep;
  int _n;
  int cnt = 0;
  _Dcomplex* inptr = 0;
  for (int p = 0; p <= nnly2; p += _step) {
   _n = n >> p;
   inptr = in + (n - _n);
   tmp1 = mem1;
   tmp2 = inptr;
   tmp3 = mem2;

   for (int step2 = 0; step2 < nly[p]; step2++) {
    step = _n >> (step2 + 1);
    tstep = step << 1;
    if (lyr[p][step2] == _n) {
     int i = 0;
     int end = lyf[p][step2] << 1;
     for (; i < end; i += tstep) {
      for (int offset = 0; offset < step; offset++) {
       _Dcomplex t = fftmath::cexpb[i << p] * tmp2[i + step + offset];
       tmp1[((i + _n) >> 1) + offset] = tmp2[i + offset] - t;
      }
     }
     for (; i < _n; i += tstep) {
      for (int offset = 0; offset < step; offset++) {
       _Dcomplex t = fftmath::cexpb[i << p] * tmp2[i + step + offset];
       tmp1[(i >> 1) + offset] = tmp2[i + offset] + t;
       tmp1[((i + _n) >> 1) + offset] = tmp2[i + offset] - t;
      }
     }
    }
    else {
     int end = lyr[p][step2] << 1;
     for (int i = (lyf[p][step2] << 1); i < end; i += tstep) {
      for (int offset = 0; offset < step; offset++) {
       _Dcomplex t = fftmath::cexpb[i << p] * tmp2[i + step + offset];
       tmp1[(i >> 1) + offset] = tmp2[i + offset] + t;
      }
     }
    }
	if (tmp2 == inptr) {
     tmp2 = tmp1;
     tmp1 = tmp3;
    }
    else {
     sw = tmp1;
     tmp1 = tmp2;
     tmp2 = sw;
    }
   }
int end = lyr[p][nly[p] - 1];
   for (int i = lyf[p][nly[p] - 1]; i < end; i++) out[cnt++] = tmp2[i];
  }
  return out;
 }
double* operator()(_Dcomplex *in) {
  out = ft(in);
  for (int b = 0; b < outsize; b++)
   que[qptr][b] = cabs(out[b]);
  for (int b = 0; b < outsize; b++) {
   dout[b] = 0.;
   int c = qptr;
   do {
    dout[b] += que[c][b];
    c += ovldata[b];
    if (c >= overlap) c -= overlap;
   } while (c != qptr);
   dout[b] *= dovldata[b];
  }
  qptr++;
  if (qptr == overlap) qptr = 0;
  return dout;
 }
virtual ~lftcalc() {
  for (int i = 0; i < nnly; i++) {
   free(lyf[i]);
   free(lyr[i]);
  }
  for (int i = 0; i < overlap; i++) free(que[i]);
  free(que);
  free(lfreq);
  free(rfreq);
  free(centerfreq);
  free(lyf);
  free(lyr);
  free(nly);
  free(mem1);
  free(mem2);
  free(out);
  free(dout);
  free(ovldata);
  free(dovldata);
 }
};
class fftcalc : virtual public fftmath {
private:
 _Dcomplex *sw, *sw2;
public:
 template<typename ...V>
 fftcalc(int overlap, int n, V...arg) : fftmath(overlap, n) {
  if (overlap < 1) throw exception("overlap must be positive integer");
  if (n&(n - 1)) throw exception("n must be power of 2");
  if (n < 2) throw exception("n must be larger than 1");
  out = (_Dcomplex*)malloc(sizeof(_Dcomplex)*n);
  outsize = n >> 1;
  dout = (double*)malloc(sizeof(double)*outsize);
  lfreq = (double*)malloc(sizeof(double)*outsize);
  rfreq = (double*)malloc(sizeof(double)*outsize);
  centerfreq = (double*)malloc(sizeof(double)*outsize);
  lfreq[0] = 0.;
  double step = 1. / (double)outsize;
  rfreq[0] = step;
  centerfreq[0] = sqrt(lfreq[0] * rfreq[0]);
  for (int i = 1; i < outsize; i++) {
   lfreq[i] = rfreq[i - 1];
   rfreq[i] = lfreq[i] + step;
   centerfreq[i] = (lfreq[i] + rfreq[i])*0.5;
  }
 }
 _Dcomplex* ft(_Dcomplex *in) {
  sw2 = out;
  for (int step = n >> 1; step > 1; step >>= 1) {
   int tstep = step << 1;
   for (int offset = 0; offset < step; offset++) {
    for (int i = 0; i < n; i += tstep) {
     _Dcomplex t = cexpb[i] * in[i + step + offset];
     sw2[(i >> 1) + offset] = in[i + offset] + t;
     sw2[((i + n) >> 1) + offset] = in[i + offset] - t;
    }
   }
   sw = in; in = sw2; sw2 = sw;
  }
  for (int i = 0; i < n; i += 2) sw2[(i >> 1)] = in[i] + cexpb[i] * in[i + 1];
  return sw2;
 }
 double *operator()(_Dcomplex *in) {
  sw2 = ft(in);
  for (int i = 0; i < outsize; i++)
   dout[i] = cabs(sw2[i]);
  return dout;
 }
 int getoutsize() {
  return outsize;
 }
 virtual ~fftcalc() {
  free(out);
  free(dout);
  free(lfreq);
  free(rfreq);
  free(centerfreq);
 }
};
class logscalefilterbank {
private:
 int nfilter, ninput;
 double **melconv;
 double* out;
 double* centerfreq;
public:
 logscalefilterbank(fftmath &input, int nfilter) :ninput(input.getoutsize()), nfilter(nfilter) {
  double* lfreq = input.getlfreq();
  double *rfreq = input.getrfreq();
  centerfreq = (double*)malloc(sizeof(double)*nfilter);
  out = (double*)malloc(sizeof(double)*nfilter);
  melconv = (double**)malloc(sizeof(double*)*nfilter);
  for (int a = 0; a < nfilter; a++) {
   melconv[a] = (double*)malloc(sizeof(double)*ninput);
   centerfreq[a] = pow(10, -4.*(double)(nfilter - 1 - a) / (double)nfilter);
   for (int b = 0; b < ninput; b++) {
    if (centerfreq[a] < lfreq[b] && centerfreq[a] < rfreq[b])
     melconv[a][b] = pow(0.0001, (lfreq[b] - centerfreq[a]) / centerfreq[a]) - pow(0.0001, (rfreq[b] - centerfreq[a]) / centerfreq[a]);
    else if (centerfreq[a] > lfreq[b] && centerfreq[a] < rfreq[b])
     melconv[a][b] = 2 - pow(0.0001, (centerfreq[a] - lfreq[b]) / centerfreq[a]) - pow(0.0001, (rfreq[b] - centerfreq[a]) / centerfreq[a]);
    else melconv[a][b] = pow(0.0001, (centerfreq[a] - rfreq[b]) / centerfreq[a]) - pow(0.0001, (centerfreq[a] - lfreq[b]) / centerfreq[a]);
   }
  }
 }
 double* operator()(double *in) {
  for (int a = 0; a < nfilter; a++) {
   out[a] = 0;
   for (int b = 0; b < ninput; b++) {
    if (melconv[a][b] > 1e-06) out[a] += in[b] * melconv[a][b];
   }
  }
  return out;
 }
 int getoutsize() {
  return nfilter;
 }
 double* getcenterfreq() {
  return centerfreq;
 }
 ~logscalefilterbank() {
  free(out);
  free(centerfreq);
 }

};

class vad {
private:
 double maxn;
 double *centerfreq, ef, em;
 int fftsize, fbsize;
 int count = 0;
 double eest, tmp, q;
 double energyupdate, energyredution, energyratio, w;
 int speechframe, hangover;
 int sp = 0, ho = 0;
 double log(double a) const noexcept {
  union { double d; long long x; } u = { a };
  return (u.x - 4606921278410026770) * 1.539095918623324e-16;
 }
public:
 double cbin(int i) {
  if (i < 0) return (centerfreq[0] + i * (centerfreq[1] - centerfreq[0]))* (double)fftsize;
  else if (i >= fbsize) return (centerfreq[fbsize - 1] + (i - fbsize + 1)*(centerfreq[fbsize - 1] - centerfreq[fbsize - 2]))* (double)fftsize;
  else return centerfreq[i] * (double)fftsize;
 }
 vad(fftmath& fft, double* centerfreq, int nfilter, double w = 1000., double energyredution = 100., double energyupdate = 20., double energyratio = 4.5, int speechframe = 3, int hangover = 6) :fbsize(nfilter), fftsize(fft.getinsize()), w(w), centerfreq(centerfreq), energyredution(energyredution), energyupdate(energyupdate), energyratio(energyratio), speechframe(speechframe), hangover(hangover) {
  maxn = 0;
  for (int i = 1; i <= fbsize; i++) maxn += (cbin(i + 1) - cbin(i - 1) + 2) *0.5;//YMAX[m,k] always 1.f
  maxn = log(maxn);
  em = 0;
 }
 vad(fftmath& fft, double w = 1000., double energyredution = 100., double energyupdate = 20., double energyratio = 4.5, int speechframe = 3, int hangover = 6) :vad(fft, fft.getcenterfreq(), fft.getoutsize(), w, energyredution, energyupdate, energyratio, speechframe, hangover) {}

 bool operator() (double* in) {
  bool ret;
  if (count < 10) {
   if (count == 0) {
    eest = 0;
    for (int i = 0; i < fbsize; i++)
     eest += in[i];
    eest = log(eest);
   }
	else {
    tmp = 0;
    for (int i = 0; i < fbsize; i++) tmp += in[i];
    eest = (eest + log(tmp))*0.5;
   }
   if (eest < (6. / 9.)*maxn) q = 32;
   else if (eest < (7. / 9.)*maxn) q = 64;
   else q = 128;
  }
  ef = 0.;
  for (int i = 0; i < fbsize; i++) ef += in[i];
  ef = q * log(1. + ef / w);
  if (!count) em = ef;
  else if (count < 10) em = (em + ef)*0.5;
  if (ef - em > energyratio) {
   sp++;
   ret = true;
}
  else {
   if (sp > speechframe) ho = hangover;
   sp = 0;
   if (ho == 0) ret = false;
   else ho--, ret = true;
  }
  if (ef - em < energyupdate) em = em + (ef - em) / energyredution;
  count++;
  return ret;
 }
};
class wavreader {
private:
 void* data;
 int smpl, nbit, count;
 _Dcomplex* dout = 0;
 int ndout = 0;
 double noiselevel = 0;
public:
 wavreader(const char* filename) {
  int hdr[11];
  FILE *f = fopen(filename, "rb");
  if (f != 0) {
   fread(hdr, 4, 11, f);
   smpl = hdr[6];
   nbit = hdr[8] >> 19;
   count = hdr[10] / nbit;
   data = reinterpret_cast<void*>(malloc(hdr[10]));
   fread(data, nbit, count, f);
   fclose(f);
  }
  else throw exception("file not exist");
 }
 void setnoiselevel(double db) {
  noiselevel = db / 32768.;
 }
 double createnoise(int n) {
  return (double)((((n + 518128927) * 1511911101 + 151211161) & 65535) - 32768)*noiselevel;
 }
 _Dcomplex* operator()(int offset, int n) {
  if (n > ndout) {
   ndout = n;
   dout = (_Dcomplex*)realloc((void*)dout, n * sizeof(_Dcomplex));
  }
  _Dcomplex *ptr1 = dout, *ptr2 = dout + n;
  switch (nbit) {
  case 1:
  {
   char* tmp = (char*)data + offset;
   char* tmp2 = (char*)data + count;
   while (ptr1 != ptr2) {
    if (tmp < tmp2) ptr1->_Val[0] = ((double)(*(tmp++) + 128)) / 128. + createnoise((int)tmp);
    else ptr1->_Val[0] = 0.;
    ptr1->_Val[1] = 0.;
    ptr1++;
   }
  }
  break;
  case 2:
  {
   short* tmp = (short*)data + offset;
   short* tmp2 = (short*)data + count;
   while (ptr1 != ptr2) {
    if (tmp < tmp2) ptr1->_Val[0] = ((double)(*(tmp++))) / 32768. + createnoise((int)tmp);
    else ptr1->_Val[0] = 0.;
    ptr1->_Val[1] = 0.;
    ptr1++;
   }
  }
  break;
  case 4:
  {
   long* tmp = (long*)data + offset;
   long* tmp2 = (long*)data + count;
   while (ptr1 != ptr2) {
    if (tmp < tmp2) ptr1->_Val[0] = ((double)(*(tmp++))) / 2147483648. + createnoise((int)tmp);
    else ptr1->_Val[0] = 0.;
    ptr1->_Val[1] = 0.;
    ptr1++;
   }
  }
  break;
  default:
   throw exception("invaild nbit");
   break;
  }
  return dout;
 }
 int getsize() {
  return count;
 }
};
class bmpwriter {
private:
 unsigned char *out = 0;
 const unsigned char red[82] = { 255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,240,225,210,195,180,165,150,135,120,105,90,75,60,45,30,15,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10,20,30,40,50,60,70,80,90,100,110,120,130,122,114,106,98,90,82,74,66,58,50,42,34,26,18,10,2,0 };
 const unsigned char grn[82] = { 255,240,225,210,195,180,165,150,135,120,105,90,75,60,45,30,15,0,15,30,45,60,75,90,105,120,135,150,165,180,195,210,225,240,255,240,225,210,195,180,165,150,135,120,105,90,75,60,45,30,15,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
 const unsigned char blu[82] = { 255,240,225,210,195,180,165,150,135,120,105,90,75,60,45,30,15,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,15,30,45,60,75,90,105,120,135,150,165,180,195,210,225,240,255,255,255,255,255,255,255,255,255,255,255,255,255,255,240,225,210,195,180,165,150,135,120,105,90,75,60,45,30,15,0 };
 const unsigned char b1[12] = { 0,0, 0,0, 54,0,0,0 ,40,0,0,0 };
 const unsigned char b2[28] = { 1,0, 24,0 };
 const unsigned char b3[3] = { 0,0,0 };
 FILE *f = 0;
 int w, h;
public:
 bmpwriter(const char* filename, int width, int height) :w(width), h(height) {
  out = (unsigned char*)malloc(3 * w*h);
  f = fopen(filename, "wb");
  if (f) {
   int filesize = 54 + 3 * w*h;
   fwrite("BM", 1, 2, f);
   fwrite(&filesize, 4, 1, f);
   fwrite(b1, 1, 12, f);
   fwrite(&w, 4, 1, f);
   fwrite(&h, 4, 1, f);
   fwrite(b2, 1, 28, f);
  }
  else throw exception("file write failed");
 }
 void operator()(double* data, int yindex) {
  for (int j = 0; j < w; j++) {
   int t = static_cast<int>(-20.*log10(1e-20 + data[j])) + 50;
   if (t < 0) t = 0;
   if (t > 81) t = 81;
   out[(j + yindex * w) * 3 + 2] = red[t]; out[(j + yindex * w) * 3 + 1] = grn[t]; out[(j + yindex * w) * 3] = blu[t];
  }
 }
 void operator()(_Dcomplex* data, int yindex) {
  for (int j = 0; j < w; j++) {
   int t = static_cast<int>(-20.*log10(1e-20 + cabs(data[j]))) + 50;
   if (t < 0) t = 0;
   if (t > 81) t = 81;
   out[(j + yindex * w) * 3 + 2] = red[t]; out[(j + yindex * w) * 3 + 1] = grn[t]; out[(j + yindex * w) * 3] = blu[t];
  }
 }
 ~bmpwriter() {
  for (int i = 0; i < h; i++) {
   fwrite(out + (w*(h - i - 1) * 3), 3, w, f);
   fwrite(b3, 1, (4 - ((w * 3) & 3)) & 3, f);
  }
  fclose(f);
  free(out);
 }
};

template<typename T>
class spectogramgenerator {
private:
 int t1 = 0, t2 = 0, t3 = 0;
public:
 template<typename ...V>
 void operator()(const char* infile, const char* outfile, int nfilter, int overlap, int n, V... arg) {
  t1 = 0; t2 = 0; t3 = 0;
  t1 = GetTickCount();
  int t21, t31, t41;
  wavreader wr(infile);
  T ft(overlap, n, arg...);
  logscalefilterbank fb(ft, nfilter);
  int shift = n / overlap;
  int size = wr.getsize();
  int y = size / shift;
  bmpwriter bm(outfile, nfilter, y);
  _Dcomplex *tmp;
  double *tmp2, *tmp3;
  for (int i = 0; i < y; i++) {
   tmp = wr(i*shift, n);
   t21 = GetTickCount();
   tmp2 = ft(tmp);
   t31 = GetTickCount();
   tmp3 = fb(tmp2);
   t41 = GetTickCount();
   bm(tmp3, i);
   t2 += t31 - t21;
   t3 += t41 - t31;
  }
  t1 = GetTickCount() - t1;
 }
 int gett1() {
  return t1;
 }
 int gett2() {
  return t2;
 }
 int gett3() {
  return t3;
 }
};

template<typename T>
class spectogramgenerator_nofb {
private:
 int t1 = 0, t2 = 0, t3 = 0;
public:
 template<typename ...V>
 void operator()(const char* infile, const char* outfile, int overlap, int n, V... arg) {
  t1 = 0; t2 = 0; t3 = 0;
	t1 = GetTickCount();
  int t21, t31, t41;
  wavreader wr(infile);
  T ft(overlap, n, arg...);
  int shift = n / overlap;
  int size = wr.getsize();
  int x = ft.getoutsize();
  int y = size / shift;
  bmpwriter bm(outfile, x, y);
  _Dcomplex *tmp;
  double *tmp2, *tmp3;
  for (int i = 0; i < y; i++) {
   tmp = wr(i*shift, n);
   t21 = GetTickCount();
   tmp2 = ft(tmp);
   t31 = GetTickCount();
   t41 = GetTickCount();
   bm(tmp2, i);
   t2 += t31 - t21;
   t3 += t41 - t31;
  }
  t1 = GetTickCount() - t1;
 }
 int gett1() {
  return t1;
 }
int gett2() {
  return t2;
 }
 int gett3() {
  return t3;
 }
};

template<typename T>
class voicedetector {
private:
 int t1 = 0, t2 = 0, t3 = 0, t4 = 0;
 bool* out;
 int y;
public:
 template<typename ...V>
 bool* operator()(const char* infile, double noiselevel, int nfilter, int overlap, int n, V... arg) {
  t1 = 0; t2 = 0; t3 = 0; t4 = 0;
  t1 = GetTickCount();
  int t21, t31, t41, t51;
  wavreader wr(infile);
  wr.setnoiselevel(noiselevel);
  T ft(overlap, n, arg...);
  logscalefilterbank fb(ft, nfilter);
  vad va(ft, fb.getcenterfreq(), nfilter);
  int shift = n / overlap;
  int size = wr.getsize();
  y = size / shift;
  out = (bool*)malloc(sizeof(bool)*y);
  _Dcomplex *tmp;
  double *tmp2, *tmp3;
  for (int i = 0; i < y; i++) {
   tmp = wr(i*shift, n);
   t21 = GetTickCount();
   tmp2 = ft(tmp);
   t31 = GetTickCount();
   tmp3 = fb(tmp2);
   t41 = GetTickCount();
   out[i] = va(tmp3);
   t51 = GetTickCount();
   t2 += t31 - t21;
   t3 += t41 - t31;
   t4 += t51 - t41;
  }
  t1 = GetTickCount() - t1;
  return out;
 }
 double compare(bool* in) {
  int diff = 0;
  for (int i = 0; i < y; i++)
   diff += (in[i] ^ out[i]);
  return (double)(y - diff) / (double)y;
 }
 int gett1() {
  return t1;
 }
 int gett2() {
  return t2;
 }
 int gett3() {
  return t3;
 }
 int gett4() {
  return t4;
 }
 ~voicedetector() {
  free(out);
 }
};
template<typename T>
class voicedetector_nofb {
private:
 int t1 = 0, t2 = 0/*, t3 = 0*/, t4 = 0;
 bool* out;
 int y;
public:
 template<typename ...V>
 bool* operator()(const char* infile, double noiselevel, int overlap, int n, V... arg) {
  t1 = 0; t2 = 0; /*t3 = 0;*/ t4 = 0;
  t1 = GetTickCount();
  int t21, t31, t41, t51;
  wavreader wr(infile);
  wr.setnoiselevel(noiselevel);
  T ft(overlap, n, arg...);
  vad va(ft, ft.getcenterfreq(), ft.getoutsize());
  int shift = n / overlap;
  int size = wr.getsize();
  y = size / shift;
  out = (bool*)malloc(sizeof(bool)*y);
  _Dcomplex *tmp;
  double *tmp2, *tmp3;
  for (int i = 0; i < y; i++) {
   tmp = wr(i*shift, n);
   t21 = GetTickCount();
   tmp2 = ft(tmp);
   t31 = GetTickCount();
   //t41 = GetTickCount();
   out[i] = va(tmp2);
   t51 = GetTickCount();
   t2 += t31 - t21;
   //t3 += t41 - t31;
   t4 += t51 - t31;
  }
  t1 = GetTickCount() - t1;
  return out;
 }
 double compare(bool* in) {
  int diff = 0;
  for (int i = 0; i < y; i++)
   diff += (in[i] ^ out[i]);
  return (double)(y - diff) / (double)y;
 }
 int gett1() {
  return t1;
 }
 int gett2() {
  return t2;
 }
 int gett3() {
  return 0;// t3;
 }
 int gett4() {
  return t4;
 }
 ~voicedetector_nofb() {
  free(out);
 }
};
void main_spectogram(int argc, char* argv[]) {
 spectogramgenerator<fftcalc> fft;
 spectogramgenerator_nofb<lftcalc> lftn;
 spectogramgenerator_nofb<fftcalc> fftn;
 char name[100];
 try {
  for (int i = 1; i < argc; i++) {
   sprintf(name, "%s_fft.bmp", argv[i]);
   fft(argv[i], name, 70, 4, 4096);
   printf("fft %d %d %d\n", fft.gett1(), fft.gett2(), fft.gett3());
   sprintf(name, "%s_lftn.bmp", argv[i]);
   lftn(argv[i], name, 4, 4096, 16, 4);
   printf("lftn %d %d %d\n", lftn.gett1(), lftn.gett2(), lftn.gett3());
   sprintf(name, "%s_fftn.bmp", argv[i]);
   fftn(argv[i], name, 4, 4096);
   printf("fftn %d %d %d\n", fftn.gett1(), fftn.gett2(), fftn.gett3());
  }
 }
 catch (exception e) {
  printf("error occured : %s\n", e.what());
 }
}
void main_vad(int argc, char* argv[]) {
 int fftt1 = 0, fftt2 = 0, fftt3 = 0, fftt4 = 0;
 int lftnt1 = 0, lftnt2 = 0, lftnt3 = 0, lftnt4 = 0;
 int fftnt1 = 0, fftnt2 = 0, fftnt3 = 0, fftnt4 = 0;
 double ffta[8] = { 0,0,0,0,0,0,0,0 }, lftna[8] = { 0,0,0,0,0,0,0,0 }, fftna[8] = { 0,0,0,0,0,0,0,0 };
 try {
  bool *t;
  for (int i = 1; i < argc; i++) {
   voicedetector<fftcalc> fft;
   voicedetector_nofb<lftcalc> lftn;
   voicedetector_nofb<fftcalc> fftn;
   double p; int ip;
   t = fft(argv[i], 0, 70, 4, 4096);
   fftt1 += fft.gett1(); fftt2 += fft.gett2(); fftt3 += fft.gett3(); fftt4 += fft.gett4();
   p = 1e-3; ip = 0;
   for (; p < 10; p *= 3.16227766017, ip++) {
    voicedetector<fftcalc> fft2;
    fft2(argv[i], p, 70, 4, 4096);
    fftt1 += fft2.gett1(); fftt2 += fft2.gett2(); fftt3 += fft2.gett3(); fftt4 += fft2.gett4();
    ffta[ip] += fft2.compare(t);
   }
   t = lftn(argv[i], 0, 4, 4096, 16, 4);
   lftnt1 += lftn.gett1(); lftnt2 += lftn.gett2(); lftnt3 += lftn.gett3(); lftnt4 += lftn.gett4();
    p = 1e-3;  ip = 0;
   for (; p < 10; p *= 3.16227766017, ip++) {
    voicedetector_nofb<lftcalc> lftn2;
    lftn2(argv[i], p, 4, 4096, 16, 4);
    lftnt1 += lftn2.gett1(); lftnt2 += lftn2.gett2(); lftnt3 += lftn2.gett3(); lftnt4 += lftn2.gett4();
    lftna[ip] += lftn2.compare(t);
   }
   t = fftn(argv[i], 0, 4, 4096);
   fftnt1 += fftn.gett1(); fftnt2 += fftn.gett2(); fftnt3 += fftn.gett3(); fftnt4 += fftn.gett4();
   p = 1e-3; ip = 0;
   for (; p < 10; p *= 3.16227766017, ip++) {
    voicedetector_nofb<fftcalc> fftn2;
    fftn2(argv[i], p, 4, 4096);
    fftnt1 += fftn2.gett1(); fftnt2 += fftn2.gett2(); fftnt3 += fftn2.gett3(); fftnt4 += fftn2.gett4();
    fftna[ip] += fftn2.compare(t);
   }
   printf("n %d %d\n", i, argc);
   printf("fft %d %d %d %d\n", fftt1, fftt2, fftt3, fftt4);
   printf("lftn %d %d %d %d\n", lftnt1, lftnt2, lftnt3, lftnt4);
   printf("fftn %d %d %d %d\n", fftnt1, fftnt2, fftnt3, fftnt4);
    p = 1e-3; ip = 0;
   for (; p < 10; p *= 3.16227766017, ip++) 
    printf("%lf fft %lf lftn %lf fftn %lf\n", p, ffta[ip] / (double)i, lftna[ip] / (double)i, fftna[ip] / (double)i);
  }
 }
 catch (exception e) {
  printf("error occured : %s\n", e.what());
 }
}
int main(int argc, char* argv[]) {
 main_vad(argc, argv);
 system("pause");
}
