#ifndef MATH_GEO_H
#define MATH_GEO_H

typedef struct {
    double* vec;//向量本体
    unsigned int w;//向量维度数
} mgeo_向量_t;

typedef struct {
    double vec[2];//向量本体
} mgeo_二维向量_t;

typedef struct {
    double vec[3];//向量本体
} mgeo_三维向量_t;

typedef struct {
    double x;
    double y;
} mgeo_二维点_t;

typedef struct {
    double x;
    double y;
    double z;
} mgeo_三维点_t;

typedef struct {
    mgeo_二维点_t* points;
    unsigned int length;
} mgeo_二维点数组_t;


mgeo_二维向量_t 点转向量_二维(mgeo_二维点_t A, mgeo_二维点_t B);
double 向量数量积_二维(mgeo_二维向量_t A, mgeo_二维向量_t B);
double 向量取模_二维(mgeo_二维向量_t A);
double 两点距离_二维(mgeo_二维点_t A, mgeo_二维点_t B);
double 三点向量夹角计算_二维(mgeo_二维点_t A, mgeo_二维点_t B, mgeo_二维点_t C);
int 最小角度点集精简_二维(mgeo_二维点数组_t *二维点数组_参数, double 最小转弯角度);
int 最小角度点集精简_二维LL(mgeo_二维点数组_t *二维点数组_参数, double 最小转弯角度);
int RDP点集精简_二维(mgeo_二维点数组_t *二维点数组_参数, double 抽除比率);
double 点到两点连线的距离(mgeo_二维点_t C, mgeo_二维点_t A, mgeo_二维点_t B);

#endif
