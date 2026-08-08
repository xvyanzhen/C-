#include "Math_Geo.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <tgmath.h>

typedef struct
{
    unsigned int *indexs;
    unsigned int length;
    unsigned int 索引容量;

} mgeo_数组索引组_t;

#define Mgo_PI 3.14159265358979323846 // M_PI在vscode中竟然提示未定义,虽然能过编译,但这里还是手动定义一个值

int 数组索引组_追加(mgeo_数组索引组_t *索引组, unsigned int 索引)
{
    // 检查内存是否足够
    if (索引组->length + 1 >= 索引组->索引容量)
    {
        unsigned int 新索引容量 = 索引组->索引容量 + (索引组->索引容量 / 2 == 0 ? 1 : 索引组->索引容量 / 2);
        unsigned int *临时指针 = realloc(索引组->indexs, 新索引容量 * sizeof(unsigned int));
        if (临时指针 == NULL)
            return -1;
        索引组->indexs = 临时指针;
        索引组->索引容量 = 新索引容量;
    }
    // 追加索引
    索引组->indexs[索引组->length] = 索引;
    索引组->length++;
    return 0;
};
int 数组索引组_释放内存(mgeo_数组索引组_t *索引组)
{
    if (索引组->indexs != NULL)
    {
        free(索引组->indexs);
        索引组->indexs = NULL;
        索引组->length = 0;
        索引组->索引容量 = 0;
    }
    return 0;
}

int 内存重分配_删除式(void **指针, unsigned int 新大小_字节)
{
    if (指针 == NULL)
        return -1;

    if (新大小_字节 == 0)
    {
        free(*指针);
        *指针 = NULL;
        return 0;
    }

    void *临时指针 = realloc(*指针, 新大小_字节);
    if (临时指针 == NULL)
        return -1;
    *指针 = 临时指针;
    return 0;
}

/**
 * @brief 初始化数组索引组
 *
 * @return mgeo_数组索引组_t 数组索引组
 */
mgeo_数组索引组_t 数组索引组初始化()
{
    mgeo_数组索引组_t 索引组 = {
        .indexs = NULL,
        .length = 0,
        .索引容量 = 0,
    };
    return 索引组;
}

/**
 * @brief 二维点转向量
 *
 * @param A 二维点A 起点
 * @param B 二维点B 终点
 * @return mgeo_向量_t 二维点A到B的向量
 */
mgeo_二维向量_t 点转向量_二维(mgeo_二维点_t A, mgeo_二维点_t B)
{
    mgeo_二维向量_t 向量 = {
        .vec = {B.x - A.x, B.y - A.y},
    };
    return 向量;
}

/**
 * @brief 二维向量数量积
 *
 * @param A 二维向量A
 * @param B 二维向量B
 * @return double 二维向量A和B的数量积
 */
double 向量数量积_二维(mgeo_二维向量_t A, mgeo_二维向量_t B)
{
    return A.vec[0] * B.vec[0] + A.vec[1] * B.vec[1];
}

/**
 * @brief 二维向量取模
 *
 * @param A 二维向量A
 * @return double 二维向量A的模
 */
double 向量取模_二维(mgeo_二维向量_t A)
{
    return sqrt(向量数量积_二维(A, A));
}

/**
 * @brief 两点距离
 *
 * @param A 二维点A 起点
 * @param B 二维点B 终点
 * @return double 两点距离
 */
double 两点距离_二维(mgeo_二维点_t A, mgeo_二维点_t B)
{
    return sqrt((B.x - A.x) * (B.x - A.x) + (B.y - A.y) * (B.y - A.y));
}

/**
 * @brief 三点折线夹角计算
 *
 * @param A 二维点A 起点
 * @param B 二维点B 中点
 * @param C 二维点C 终点
 * @return double 三点折线夹角
 */
double 三点向量夹角计算_二维(mgeo_二维点_t A, mgeo_二维点_t B, mgeo_二维点_t C)
{
    mgeo_二维向量_t AB = 点转向量_二维(A, B);
    mgeo_二维向量_t BC = 点转向量_二维(B, C);
    double AB与BC的数量积 = 向量数量积_二维(AB, BC);
    double AB的模 = 向量取模_二维(AB);
    if (AB的模 <= 0)
        return NAN;
    double BC的模 = 向量取模_二维(BC);
    if (BC的模 <= 0)
        return NAN;

    double 余弦值 = AB与BC的数量积 / (AB的模 * BC的模);
    // 钳位到[-1,1]防止浮点误差导致acos参数越界
    if (余弦值 > 1.0)
        余弦值 = 1.0;
    if (余弦值 < -1.0)
        余弦值 = -1.0;

    return acos(余弦值);
}

/**
 * @brief C点到AB两点连线的距离
 *
 * @param C 二维点C 被计算点
 * @param A 二维点A 起点
 * @param B 二维点B 终点
 * @return double 点到两点连线的距离
 */
double 点到两点连线的距离(mgeo_二维点_t C, mgeo_二维点_t A, mgeo_二维点_t B)
{
    if (A.x == B.x && A.y == B.y)
    {
        return NAN;
    };
    mgeo_二维向量_t AB = 点转向量_二维(A, B);
    double AB的模 = 向量取模_二维(AB);
    if (AB的模 <= 0)
    {
        return NAN;
    }
    mgeo_二维向量_t AC = 点转向量_二维(A, C);
    double AC的模 = 向量取模_二维(AC);
    double AC投影长 = 向量数量积_二维(AB, AC) / AB的模;
    double 垂线长度 = sqrt(AC的模 * AC的模 - AC投影长 * AC投影长);
    return 垂线长度;
}

/**
 * @brief 最小角度点集精简,每次删除夹角小于最小转弯角度的点,并保护被删除点的两端的点不被删除,直到所有点都被删除
 *
 * @param 二维点数组 传入的二维点数组,他的内部数组会被修改
 * @param 最小转弯角度 转弯截止角度,单位为度
 * @return mgeo_二维点数组_t 精简后的二维点数组
 */
int 最小角度点集精简_二维(mgeo_二维点数组_t *二维点数组_参数, double 最小转弯角度)
{
    mgeo_二维点数组_t 二维点数组 = *二维点数组_参数;
    if (二维点数组.length == 0)
        return 0;
    if (二维点数组.length == 1)
        return 0;
    if (二维点数组.length == 2)
        return 0;

    double 最小转弯角度_弧度 = 最小转弯角度 * Mgo_PI / 180;
    // 由于起点和重点定然跳过,且会保护被删除点的两端不被删除,就是间隔删除,所以不需要全部初始化全部长度
    unsigned int 最大删除索引组长度 = 二维点数组.length / 2 + 1; // 最坏情况就是连续每间隔一个点就会有一个点被删除,这会产生 总数/2 个删除点索引
    mgeo_数组索引组_t 删除点索引组 = 数组索引组初始化();
    删除点索引组.索引容量 = 最大删除索引组长度 < 500 ? 最大删除索引组长度 : 500;
    删除点索引组.indexs = malloc(删除点索引组.索引容量 * sizeof(unsigned int));
    if (删除点索引组.indexs == NULL)
    {
        goto 错误处理;
    }

    // 初次标记删除点索引
    for (unsigned int i = 1; i < 二维点数组.length - 1; i++)
    {
        double 夹角_弧度 = 三点向量夹角计算_二维(二维点数组.points[i - 1], 二维点数组.points[i], 二维点数组.points[i + 1]);
        if (夹角_弧度 <= 最小转弯角度_弧度 || isnan(夹角_弧度))
        {
            // 检查删除点索引组是否够大
            if (删除点索引组.length >= 删除点索引组.索引容量)
            {
                unsigned int 新大小 = 删除点索引组.索引容量 * 2 < 最大删除索引组长度 ? 删除点索引组.索引容量 * 2 : 最大删除索引组长度;
                if (内存重分配_删除式((void **)&删除点索引组.indexs, 新大小 * sizeof(unsigned int)) != 0)
                {
                    goto 错误处理;
                }
                删除点索引组.索引容量 = 新大小;
            }
            删除点索引组.indexs[删除点索引组.length] = i; // 将新点附加到数组末尾
            删除点索引组.length++;                        // 删除点索引组长度增加1

            i++; // 防止删除删除点的两端的点
        }
    }

    while (删除点索引组.length > 0)
    {

        // 删除点索引组中所有点
        unsigned int 删除后二维点数组_length = 二维点数组.length - 删除点索引组.length;
        // 将点复制到缓冲区
        unsigned int 写入索引 = 0;
        for (unsigned int i = 0, j = 0; i < 二维点数组.length; i++)
        {
            if (j < 删除点索引组.length && 删除点索引组.indexs[j] == i)
            {
                j++;
                continue;
            }
            else
            {
                二维点数组.points[写入索引] = 二维点数组.points[i];
                写入索引++;
            }
        }

        // 大小更新,实际内存的长度更新会在函数执行末尾进行
        二维点数组.length = 删除后二维点数组_length;

        // 删除点标记重置
        删除点索引组.length = 0;

        // 标记删除点索引
        for (unsigned int i = 1; i < 二维点数组.length - 1; i++)
        {
            double 夹角_弧度 = 三点向量夹角计算_二维(二维点数组.points[i - 1], 二维点数组.points[i], 二维点数组.points[i + 1]);
            if (夹角_弧度 <= 最小转弯角度_弧度 || isnan(夹角_弧度))
            {
                // 检查删除点索引组是否够大
                if (删除点索引组.length >= 删除点索引组.索引容量)
                {
                    unsigned int 新大小 = 删除点索引组.索引容量 * 2 < 最大删除索引组长度 ? 删除点索引组.索引容量 * 2 : 最大删除索引组长度;
                    if (内存重分配_删除式((void **)&删除点索引组.indexs, 新大小 * sizeof(unsigned int)) != 0)
                    {
                        goto 错误处理;
                    }
                    删除点索引组.索引容量 = 新大小;
                }
                删除点索引组.indexs[删除点索引组.length] = i; // 将新点附加到数组末尾
                删除点索引组.length++;                        // 删除点索引组长度增加1

                i++; // 防止删除删除点的两端的点
            }
        }
    };

    if (内存重分配_删除式((void **)&二维点数组.points, 二维点数组.length * sizeof(mgeo_二维点_t)) != 0)
    {
        // 若果分配失败,意味着结果占用着过大的内存,不过依然可以使用,只是内存占用会增加
    }
    二维点数组_参数->points = 二维点数组.points;
    二维点数组_参数->length = 二维点数组.length;

    free(删除点索引组.indexs);

    return 0;

错误处理:
    free(删除点索引组.indexs);
    return -1;
}

/**
 * @brief 最小角度点集精简的超长点集输入版本,每次删除夹角小于最小转弯角度的点,并保护被删除点的两端的点不被删除,直到所有点都被删除.只会在点集长度500万以上,且可删除点稀疏(最小转弯角度较小)的情况下速度才可能会超过普通版本
 *
 * @param 二维点数组 传入的二维点数组,他的内部数组会被修改
 * @param 最小转弯角度 转弯截止角度,单位为度
 * @return mgeo_二维点数组_t 精简后的二维点数组
 */
int 最小角度点集精简_二维LL(mgeo_二维点数组_t *二维点数组_参数, double 最小转弯角度)
{
    mgeo_二维点数组_t 二维点数组 = *二维点数组_参数;
    if (二维点数组.length == 0)
        return 0;
    if (二维点数组.length == 1)
        return 0;
    if (二维点数组.length == 2)
        return 0;

    double 最小转弯角度_弧度 = 最小转弯角度 * Mgo_PI / 180;
    // 由于起点和终点定然跳过,且会保护被删除点的两端不被删除,就是间隔删除,所以不需要全部初始化全部长度
    unsigned int 最大删除索引组长度 = 二维点数组.length / 2 + 1;       // 最坏情况就是连续每间隔一个点就会有一个点被删除,这会产生 总数/2 个删除点索引
    unsigned int 最大检查点索引组长度 = 二维点数组.length * 2 / 3 + 1; // 最坏情况就是每三个点的中间点都被检查,这会产生 总数*2/3 个检查点索引
    mgeo_数组索引组_t 删除点索引组 = 数组索引组初始化();
    mgeo_数组索引组_t 检查点索引组 = 数组索引组初始化();

    
    删除点索引组.索引容量 = 最大删除索引组长度 < 500 ? 最大删除索引组长度 : 500;
    删除点索引组.indexs = malloc(删除点索引组.索引容量 * sizeof(unsigned int));
    if (删除点索引组.indexs == NULL)
        goto 错误处理;

    检查点索引组.索引容量 = 最大检查点索引组长度 < 500 ? 最大检查点索引组长度 : 500;
    检查点索引组.indexs = malloc(检查点索引组.索引容量 * sizeof(unsigned int));
    if (检查点索引组.indexs == NULL)
        goto 错误处理;
    // 初次标记删除点索引
    for (unsigned int i = 1; i < 二维点数组.length - 1; i++)
    {
        double 夹角_弧度 = 三点向量夹角计算_二维(二维点数组.points[i - 1], 二维点数组.points[i], 二维点数组.points[i + 1]);
        if (夹角_弧度 <= 最小转弯角度_弧度 || isnan(夹角_弧度))
        {
            // 检查删除点索引组是否够大
            if (删除点索引组.length >= 删除点索引组.索引容量)
            {
                unsigned int 新大小 = 删除点索引组.索引容量 * 2 < 最大删除索引组长度 ? 删除点索引组.索引容量 * 2 : 最大删除索引组长度;
                if (内存重分配_删除式((void **)&删除点索引组.indexs, 新大小 * sizeof(unsigned int)) != 0)
                {
                    goto 错误处理;
                }
                删除点索引组.索引容量 = 新大小;
            }
            删除点索引组.indexs[删除点索引组.length] = i; // 将新点附加到数组末尾
            删除点索引组.length++;                        // 删除点索引组长度增加1

            i++; // 防止删除删除点的两端的点
        }
    }

    while (删除点索引组.length > 0)
    {
        检查点索引组.length = 0;
        // 删除点索引组中所有点,并标记需要检查的点
        int 删除后长度 = 二维点数组.length - 删除点索引组.length;
        if (删除后长度 <= 0)
            goto 错误处理;
        unsigned int 写入索引 = 0;
        // 将点复制到缓冲区
        for (unsigned int i = 0, j = 0; i < 二维点数组.length; i++)
        {
            if (j < 删除点索引组.length && 删除点索引组.indexs[j] == i)
            {
                j++;
                // 标记需要检查的点

                // 添加前驱点
                unsigned int 前驱点索引 = 写入索引 - 1;
                // 检查内存是否充足,且前驱点索引不是 当前点索引减1 则添加前驱点,用来防止重复添加
                if (前驱点索引 > 0 && (检查点索引组.length == 0 || 检查点索引组.indexs[检查点索引组.length - 1] != 写入索引 - 1))
                {
                    if (检查点索引组.length >= 检查点索引组.索引容量)
                    {
                        unsigned int 新大小 = 检查点索引组.索引容量 * 2 < 最大检查点索引组长度 ? 检查点索引组.索引容量 * 2 : 最大检查点索引组长度;
                        if (内存重分配_删除式((void **)&检查点索引组.indexs, 新大小 * sizeof(unsigned int)) != 0)
                        {
                            goto 错误处理;
                        }
                        检查点索引组.索引容量 = 新大小;
                    }
                    检查点索引组.indexs[检查点索引组.length] = 前驱点索引; // 将新点在新数组对应的索引附加到数组末尾
                    检查点索引组.length++;                                 // 检查点索引组长度增加1
                }

                // 添加后继点
                unsigned int 后继点索引 = 写入索引;
                if (后继点索引 != (unsigned int)(删除后长度-1)) // 防止添加最后一个点的索引值
                {
                    { // 检查内存是否充足
                        if (检查点索引组.length >= 检查点索引组.索引容量)
                        {
                            unsigned int 新大小 = 检查点索引组.索引容量 * 2 < 最大检查点索引组长度 ? 检查点索引组.索引容量 * 2 : 最大检查点索引组长度;
                            if (内存重分配_删除式((void **)&检查点索引组.indexs, 新大小 * sizeof(unsigned int)) != 0)
                            {
                                goto 错误处理;
                            }
                            检查点索引组.索引容量 = 新大小;
                        }
                        检查点索引组.indexs[检查点索引组.length] = 写入索引; // 将新点在新数组对应的索引附加到数组末尾
                        检查点索引组.length++;
                    } // 检查点索引组长度增加1
                }
                continue;
            }
            else
            {
                二维点数组.points[写入索引] = 二维点数组.points[i];
                写入索引++;
            }
        }

        // 大小更新,实际内存的长度更新会在函数执行末尾进行
        二维点数组.length = (unsigned int)删除后长度;
        // 删除点标记重置
        删除点索引组.length = 0;

        // 检测要删除的点
        unsigned int 下一个跳过位置 = 0;
        for (int i = 0; i < 检查点索引组.length; i++)
        {
            if (检查点索引组.indexs[i] == 下一个跳过位置)
            {
                // 由于所有i中满足 i === 下一个跳过位置 的i只会有一个,所以这里直接跳过即可
                continue;
            }
            double 夹角_弧度 = 三点向量夹角计算_二维(二维点数组.points[检查点索引组.indexs[i] - 1], 二维点数组.points[检查点索引组.indexs[i]], 二维点数组.points[检查点索引组.indexs[i] + 1]);
            if (夹角_弧度 <= 最小转弯角度_弧度 || isnan(夹角_弧度))
            {
                // 检查删除点索引组是否够大
                if (删除点索引组.length >= 删除点索引组.索引容量)
                {
                    unsigned int 新大小 = 删除点索引组.索引容量 * 2 < 最大删除索引组长度 ? 删除点索引组.索引容量 * 2 : 最大删除索引组长度;
                    if (内存重分配_删除式((void **)&删除点索引组.indexs, 新大小 * sizeof(unsigned int)) != 0)
                    {
                        goto 错误处理;
                    }
                    删除点索引组.索引容量 = 新大小;
                }
                删除点索引组.indexs[删除点索引组.length] = 检查点索引组.indexs[i]; // 将新点附加到数组末尾
                删除点索引组.length++;

                下一个跳过位置 = 检查点索引组.indexs[i] + 1;
            }
        }

    };

    if (内存重分配_删除式((void **)&二维点数组.points, 二维点数组.length * sizeof(mgeo_二维点_t)) != 0)
    {
        // 若果分配失败,意味着结果占用着过大的内存,不过依然可以使用,只是内存占用会增加
    }
    // 参数同步
    二维点数组_参数->length = 二维点数组.length;
    二维点数组_参数->points = 二维点数组.points;

    free(检查点索引组.indexs);
    free(删除点索引组.indexs);

    return 0;

错误处理:
    free(检查点索引组.indexs);
    free(删除点索引组.indexs);
    return -1;
}

/**
 * @brief RDP点集精简
 *
 * @param 二维点数组_参数 二维点数组
 * @param 抽除比率 抽除比率
 * @return int 0 成功 -1 失败
 */
int RDP点集精简_二维(mgeo_二维点数组_t *二维点数组_参数, double 抽除比率)
{
    mgeo_二维点数组_t 二维点数组 = *二维点数组_参数;
    if (二维点数组.points == NULL)
        return -1;
    if (二维点数组.length == 0)
        return 0;
    if (二维点数组.length == 1)
        return 0;
    if (二维点数组.length == 2)
        return 0;

    double *中间点距离 = malloc((二维点数组.length - 2) * sizeof(double));
    mgeo_数组索引组_t 循环范围栈 = 数组索引组初始化();
    mgeo_数组索引组_t 保留点索引栈 = 数组索引组初始化();

    // 初始化存储中间点距离的数组
    if (中间点距离 == NULL)
    {
        goto 错误处理;
    }
    unsigned int 中间点距离_索引容量 = 二维点数组.length - 2;

    // 初始化循环范围栈
    {
        unsigned int 初始分配的索引容量 = 二维点数组.length / 2 < 500 ? 二维点数组.length / 2 : 500;
        循环范围栈.indexs = malloc(初始分配的索引容量 * sizeof(unsigned int));
        if (循环范围栈.indexs == NULL)
        {
            goto 错误处理;
        }
        循环范围栈.索引容量 = 初始分配的索引容量;
    }

    // 初始化保留点索引栈
    {
        unsigned int 初始分配的索引容量 = 二维点数组.length / 2 < 500 ? 二维点数组.length / 2 : 500;
        保留点索引栈.indexs = malloc(初始分配的索引容量 * sizeof(unsigned int));
        if (保留点索引栈.indexs == NULL)
        {
            goto 错误处理;
        }
        保留点索引栈.索引容量 = 初始分配的索引容量;
    }
    if (数组索引组_追加(&保留点索引栈, 0) != 0)
        goto 错误处理;
    if (数组索引组_追加(&循环范围栈, 0) != 0)
        goto 错误处理;
    if (数组索引组_追加(&循环范围栈, 二维点数组.length - 1) != 0)
        goto 错误处理;

    while (循环范围栈.length > 0)
    {

        unsigned int 中间点距离_已用索引 = 0;
        // 取出循环范围
        unsigned int 左边界 = 循环范围栈.indexs[循环范围栈.length - 2];
        unsigned int 右边界 = 循环范围栈.indexs[循环范围栈.length - 1];
        循环范围栈.length -= 2;
        if (右边界 - 左边界 == 1)
        {
            if (数组索引组_追加(&保留点索引栈, 右边界) != 0)
                goto 错误处理; // 保留右边界 ,最初的左边界已添加,由于在后续循环中本次循环的右边界就是后续循环的左边界,所以不添加左边界
            continue;
        }
        mgeo_二维向量_t 起止点向量 = 点转向量_二维(二维点数组.points[左边界], 二维点数组.points[右边界]);
        double 起止点距离 = 向量取模_二维(起止点向量);
        mgeo_二维点_t 起点 = 二维点数组.points[左边界];
        if (起止点距离 <= 0) // 这代表着起止点重合,无法判断是否需要保留,以距离起点最远的点为分界线再次分割为两个范围
        {
            for (unsigned int i = 左边界 + 1; i < 右边界; i++)
            {
                中间点距离[中间点距离_已用索引] = 两点距离_二维(起点, 二维点数组.points[i]);
                中间点距离_已用索引++;
            }
            unsigned int 最大距离索引 = 0;
            double 最大值 = 中间点距离[0];
            for (unsigned int i = 1; i < 中间点距离_已用索引; i++)
            {
                if (中间点距离[i] > 最大值)
                {
                    最大值 = 中间点距离[i];
                    最大距离索引 = i;
                }
            }
            // 栈会后进先出,先压入后半段范围,再压入前半段范围,可以保证下次先扫描前半段的范围,保证添加的保留点索引从小到大排序
            // 后半段范围压入
            if (数组索引组_追加(&循环范围栈, 左边界 + 最大距离索引 + 1) != 0)
                goto 错误处理;
            if (数组索引组_追加(&循环范围栈, 右边界) != 0)
                goto 错误处理;
            // 前半段范围压入
            if (数组索引组_追加(&循环范围栈, 左边界) != 0)
                goto 错误处理;
            if (数组索引组_追加(&循环范围栈, 左边界 + 最大距离索引 + 1) != 0)
                goto 错误处理;
            continue;
        }
        for (unsigned int i = 左边界 + 1; i < 右边界; i++)
        {
            mgeo_二维向量_t AC = 点转向量_二维(起点, 二维点数组.points[i]);
            double AC的模 = 向量取模_二维(AC);
            double AC投影长 = 向量数量积_二维(起止点向量, AC) / 起止点距离;
            double 垂线长度 = sqrt(AC的模 * AC的模 - AC投影长 * AC投影长);
            中间点距离[中间点距离_已用索引] = 垂线长度;
            中间点距离_已用索引++;
        }
        unsigned int 最大距离索引 = 0;
        double 最大值 = 中间点距离[0];
        for (unsigned int i = 1; i < 中间点距离_已用索引; i++)
        {
            if (中间点距离[i] > 最大值)
            {
                最大值 = 中间点距离[i];
                最大距离索引 = i;
            }
        }
        if (最大值 > 抽除比率 * 起止点距离)
        {
            // 栈会后进先出,先压入后半段范围,再压入前半段范围,可以保证下次先扫描前半段的范围,保证添加的保留点索引从小到大排序
            // 后半段范围压入
            if (数组索引组_追加(&循环范围栈, 左边界 + 最大距离索引 + 1) != 0)
                goto 错误处理;
            if (数组索引组_追加(&循环范围栈, 右边界) != 0)
                goto 错误处理;

            // 前半段范围压入
            if (数组索引组_追加(&循环范围栈, 左边界) != 0)
                goto 错误处理;
            if (数组索引组_追加(&循环范围栈, 左边界 + 最大距离索引 + 1) != 0)
                goto 错误处理;
        }
        else
        {
            // 保留右边界
            if (数组索引组_追加(&保留点索引栈, 右边界) != 0)
                goto 错误处理;
        }
    }
    free(中间点距离);
    中间点距离 = NULL;

    for (unsigned int i = 0; i < 保留点索引栈.length; i++)
    {
        二维点数组.points[i] = 二维点数组.points[保留点索引栈.indexs[i]];
    }

    if (内存重分配_删除式((void **)&二维点数组.points, 保留点索引栈.length * sizeof(mgeo_二维点_t)) != 0)
    {
        // 若果分配失败,意味着结果占用着过大的内存,不过依然可以使用,只是内存占用会增加
    }
    二维点数组_参数->points = 二维点数组.points;
    二维点数组_参数->length = 保留点索引栈.length;

    free(保留点索引栈.indexs);
    free(循环范围栈.indexs);
    return 0;

错误处理:
    free(中间点距离);
    free(保留点索引栈.indexs);
    free(循环范围栈.indexs);
    return -1;
}
