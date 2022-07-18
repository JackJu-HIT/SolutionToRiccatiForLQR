/**
 * @Function:Solution for Riccati 
 * @Date:2022-07-18 14:00:00
 * @Author:juchunyu
 * @Last modified:juchunyu
 */
#include <stdio.h>
#include <math.h>

#define SIZE  3
#define SIZE_ 1

int main(){

    double A[SIZE][SIZE]                        = {0,1,0,0,0,1,-7,-41,6};
    double B[SIZE]                              = {0,0,1};
    double C[SIZE]                              = {6,0,0};
    double Q_C[SIZE][SIZE]                      = {1,0,1,0,1,0,0,0,1};
    double R_C[SIZE_ ][SIZE_ ]                  = {2};
    double AT[SIZE][SIZE]                       = {0}; 
    double BT[SIZE]                             = {0};    
    double P_Result[SIZE][SIZE]                 = {0};
    double P[SIZE][SIZE]                        = {1,0,1,0,1,0,0,0,1};              //P的初始值设定为Q矩阵

    int    times                                = 50; //迭代最大次数限制
    double tolerance                            = 0.01; //迭代精度设置
    int    IsResult                             = 0;

    double  K[SIZE]                             = {0};

    //计算A的转置
    for(int i = 0;i < SIZE;i++){
        for(int j = 0;j < SIZE;j++){
            AT[j][i] = A[i][j];
        }
    }

    //计算B的转置
     for(int i = 0;i < SIZE;i++){
         BT[i] = B[i];
    }

    /*P_next = A'.*P.*A-A'.*P.*B.*(R+B'.*P.*B)^(-1).*B'.*P.*A+Q*/
    for(int t = 0;t < times;t++){
      
        double P_next[SIZE][SIZE]                   = {0};
        double ATP[SIZE][SIZE]                      = {0};                
        double ATPA[SIZE][SIZE]                     = {0};
        double ATPB[SIZE]                           = {0};
        double BTP[SIZE]                            = {0};
        double BTPA[SIZE]                           = {0};
        double BTPB                                 =  0;
        double RBTPB[SIZE_][SIZE_]                  = {0};
        double RBTPBINVERSE[SIZE_][SIZE_]           = {0};
        double ATPBRBTPBINVERSE[SIZE]               = {0};
        double ATPBRBTPBINVERSEBTPA[SIZE][SIZE]     = {0};
      

        printf("矩阵P的结果：\n");
        for(int i = 0;i < SIZE;++i)
        {
            for(int j = 0;j < SIZE;++j)
            {
                printf("%2.4f    ",P[i][j]);
            }
            printf("\n");
        }
        
        //计算A'*P
        for(int i = 0;i < SIZE;i++){
            for(int j = 0;j < SIZE;j++){
                for(int k = 0;k < SIZE;k++){
                    ATP[i][j] += AT[i][k] * P[k][j];
                }
            }
        }
        printf("矩阵AT的结果：\n");
        for(int i = 0;i < SIZE;++i)
        {
            for(int j = 0;j < SIZE;++j)
            {
                printf("%2.4f    ",AT[i][j]);
            }
            printf("\n");
        }
        //计算A'*P*A
        for(int i = 0;i < SIZE;i++){
            for(int j = 0;j < SIZE;j++){
                for(int k = 0;k < SIZE;k++){
                    ATPA[i][j] += ATP[i][k] * A[k][j];
                }
            }
        }

        printf("矩阵A'PA的结果：\n");
        for(int i = 0;i < SIZE;++i)
        {
            for(int j = 0;j < SIZE;++j)
            {
                printf("%2.4f    ",ATPA[i][j]);
            }
            printf("\n");
        }

        //计算A'*P*B
        for(int i = 0;i < SIZE;i++){
                for(int k = 0;k < SIZE;k++){
                    ATPB[i] += ATP[i][k] * B[k];
                }
        }
        
        printf("矩阵A'PB的结果：\n");
        for(int i = 0;i < SIZE;++i)
        {
  
            printf("%2.4f    ",ATPB[i]);
            printf("\n");
        }

        //计算B'*P
        for(int i = 0;i < SIZE;i++){
                for(int k = 0;k < SIZE;k++){
                    BTP[i] += B[k] * P[k][i];
                }
        }
        
        //计算B'*P*A
        for(int i = 0;i < SIZE;i++){
                for(int k = 0;k < SIZE;k++){
                    BTPA[i] += BTP[k] * A[k][i];
                }
        }

        printf("矩阵B'PA的结果：\n");
        for(int i = 0;i < SIZE;++i)
        {
            printf("%2.4f    ",BTPA[i]);
            printf("\n");
        }

        //计算B'*P*B
        for(int i = 0; i < SIZE; i++){
           BTPB += BTP[i] * B[i];  
        }

        printf("矩阵B'PB的结果：\n");
        printf("%2.4f    ",BTPB);
        printf("\n");

        //计算R+B'*P*B
        for(int i = 0;i < SIZE_;i++){
            for(int j = 0;j < SIZE_;j++){
                RBTPB[i][j] += R_C[i][j] + BTPB; 
            }
        }
        printf("矩阵R+B'*P*B的结果：\n");
        for(int i = 0;i < SIZE_;++i)
        {
            for(int j = 0;j < SIZE_;++j)
            {
                printf("%2.4f    ",RBTPB[i][j]);
            }
            printf("\n");
        }
      
        /******************************************开始计算(R+B'*P*B)^(-1)********************************************************/
        double Q[SIZE_][SIZE_]              = {0};
        double R[SIZE_][SIZE_]              = {0};
        double v[SIZE_]                     = {0};
        double QT[SIZE_][SIZE_]             = {0};
        double RInverse[SIZE_][SIZE_]       = {0};
        double m_in[SIZE_][SIZE_]           = {0};
        double m_inverse_[SIZE_][SIZE_]     = {0};
        double m_inverse[SIZE_][SIZE_]      = {0}; 
        int    length                       = SIZE_;  
        /**进行QR分解**/
        for(int k = 0;k < length;k++){

            /*R(1:k-1,k) = Q(:,1:k-1)’ * A(:,k)*/
            if(k >= 1){
                for(int i = 0;i < k;i++){
                    for(int j = 0;j < length;j++){
                        R[i][k] += Q[j][i]*RBTPB[j][k];
                    }
                }
            }
            /*v = A(:,k) - Q(:,1:k-1) * R(1:k-1,k)*/
            for(int j = 0;j<length;j++){
                if(k < 1){
                v[j] = RBTPB[j][0];
                } else {
                        v[j] = RBTPB[j][k];
                        for(int g = 0;g < k;g++){
                            v[j] -= R[g][k]*Q[j][g];
                        }
                }
            }

            /*R(k,k) = norm(v)*/
            for(int i = 0;i < length;i++){
                R[k][k] += v[i]*v[i];
            }

            R[k][k] = sqrt(R[k][k]);

            /*Q(:,k) = v / R(k,k)*/
            for(int i = 0;i < length;i++){
                Q[i][k] = v[i]/R[k][k];
            }  
        }
    
        //求解Q的转置
        for(int i = 0;i < length;i++){
        for(int j = 0;j<length;j++){
            QT[j][i] = Q[i][j];
        }
        }
    
        /**求解R的逆矩阵***/
        //转置
        for(int i = 0;i < length;i++){
            for(int j = 0;j < length;j++){
                m_in[j][i] = R[i][j];
            }
        }
        
    
        for(int j = 0;j < length;j++){
            //求解对角线
            m_inverse_[j][j] = 1/m_in[j][j];
            for(int i = j+1;i < length;i++){
                    double temp = 0;

                    for(int k = j;k <= i-1;k++){
                        temp += m_in[i][k]*m_inverse_[k][j];
                    }

                    m_inverse_[i][j] = -temp/m_in[i][i];
            }
        }
    
        //tranpose
        for(int i = 0;i < length;i++){
            for(int j = 0;j < length;j++){
                m_inverse[j][i] = m_inverse_[i][j];  //逆矩阵
            }
        }
    
        /*A = QR => A^(-1) = R^(-1)*Q^T*/

        for(int i = 0;i < length;i++){
            for(int j = 0;j < length;j++){
                for(int k = 0;k < length;k++){
                    RBTPBINVERSE[i][j] += m_inverse[i][k]*QT[k][j]; 
                }
            }
        }

        //打印
        printf("矩阵(R+B'*P*B)^(-1)的结果：\n");
        for(int i = 0;i < SIZE_;i++){
            for(int j = 0;j < SIZE_;++j){
                printf("%2.4f    ",RBTPBINVERSE[i][j]);
            }
            printf("\n");
        }

        /****************************************************************结束(R+B'*P*B)^(-1)求解********************************************/
        
        //计算A'*P*B*(R+B'*P*B)^(-1)
        for(int i = 0;i < SIZE;i++){
            for(int k = 0;k < SIZE_;k++){
                ATPBRBTPBINVERSE[i] += ATPB[i] * RBTPBINVERSE[k][0];
            }
        }
        printf("矩阵A'*P*B*(R+B'*P*B)^(-1)的结果：\n");
        for(int i = 0;i < SIZE;i++){
            printf("%2.4f    ",ATPBRBTPBINVERSE[i]);
            printf("\n");
        }

        //计算A'*P*B*(R+B'*P*B)^(-1)*B'*P*A
        for(int i = 0;i < SIZE;i++){
            for (int j = 0; j < SIZE; j++){
                    ATPBRBTPBINVERSEBTPA[i][j] = ATPBRBTPBINVERSE[i] *  BTPA[j];
            }
            
        }
        
        printf("矩阵A'*P*B*(R+B'*P*B)^(-1)*B'*P*A的结果：\n");
        
        for(int i = 0;i < SIZE;i++){
            for(int j = 0;j < SIZE;++j){
                printf("%2.4f    ", ATPBRBTPBINVERSEBTPA[i][j]);
            }
            printf("\n");
        }


        /*P_next = A'.*P.*A-A'.*P.*B.*(R+B'.*P.*B)^(-1).*B'.*P.*A+Q*/
        for(int i = 0;i< SIZE;i++){
            for (int j = 0; j < SIZE; j++){
               P_next[i][j] = ATPA[i][j] -  ATPBRBTPBINVERSEBTPA[i][j] + Q_C[i][j];
            }
        }

        printf("矩阵P_next的结果：\n");
        for(int i = 0;i < SIZE;i++){
            for(int j = 0;j < SIZE;++j){
                printf("%2.4f    ", P_next[i][j]);
            }
            printf("\n");
        }

        //判断是否收敛
        /*
        double P_maxCoeff = 0;
        double P_max[SIZE][SIZE] = {0};

        for(int i = 0;i < SIZE;i++){
            for(int j = 0;j < SIZE;j++){
                P_max[i][j] = (P[i][j] - P_next[i][j]) > 0 ?  (P[i][j] - P_next[i][j]):-(P[i][j] - P_next[i][j]) ;
            }
        }

        for(int i = 0;i < SIZE;i++){
            for(int j = 0;j < SIZE;j++){
                if(P_max[i][j]>P_maxCoeff)
                    P_maxCoeff = P[i][j];
            }
        }
        */

        /*
        for(int i = 0;i < SIZE;i++){
            for(int j = 0;j < SIZE;j++){
                if(P[i][j]>P_maxCoeff)
                P_maxCoeff = P[i][j];
            }
        }
        */
        /*
        for(int i = 0;i < SIZE;i++){
            for(int j = 0;j < SIZE;j++){
                P_maxCoeff += P[i][j] * P[i][j];
            }
        }
        P_maxCoeff = sqrt(P_maxCoeff);
        */
        
        double P_next_maxCoeff = 0;
       /*
        for(int i = 0;i < SIZE;i++){
            for(int j = 0;j < SIZE;j++){
                if(P_next[i][j]>P_maxCoeff)
                P_next_maxCoeff = P[i][j];
            }
        }
        */
        /*
        for(int i = 0;i < SIZE;i++){
            for(int j = 0;j < SIZE;j++){
                P_next_maxCoeff = P_next[i][j] * P_next[i][j];
            }
        }
         P_next_maxCoeff = sqrt(P_next_maxCoeff);
        */
       /*
        if((P_next_maxCoeff - P_maxCoeff) > 0){
            if((P_next_maxCoeff - P_maxCoeff) < tolerance){
                //将P_Next赋值给P_Result
                for(int i = 0;i < SIZE;i++){
                    for(int j = 0;j < SIZE;j++){
                        P_Result[i][j] = P_next[i][j];
                    }
                }
                IsResult = 1;
                break;
            }
        } 
        
        if((P_next_maxCoeff - P_maxCoeff) < 0){
            if(-(P_next_maxCoeff - P_maxCoeff) < tolerance){
                //将P_Next赋值给P_Result
                for(int i = 0;i < SIZE;i++){
                    for(int j = 0;j < SIZE;j++){
                        P_Result[i][j] = P_next[i][j];
                    }
                }
                IsResult = 1;
                break;
            }
        }
         */
        /*
        if(P_maxCoeff < tolerance){
                //将P_Next赋值给P_Result
                for(int i = 0;i < SIZE;i++){
                    for(int j = 0;j < SIZE;j++){
                        P_Result[i][j] = P_next[i][j];
                    }
                }
                IsResult = 1;
                break;
        }
        */
        //将P_Next赋值给P_Result
        for(int i = 0;i < SIZE;i++){
            for(int j = 0;j < SIZE;j++){
                P_Result[i][j] = P_next[i][j];
            }
        }
        IsResult = 1;
        //将P_Next赋值给P
        for(int i = 0;i < SIZE;i++){
            for(int j = 0;j < SIZE;j++){
                P[i][j] = P_next[i][j];
            }
        }

         /***K = RBTPBINVERSE[k][0]**/
        for(int i = 0;i < SIZE;i++){
            K[i] =  RBTPBINVERSE[0][0] * BTPA[i];
        }
    } 
    if(IsResult == 1){
            printf("矩阵P_Result的结果：\n");
            for(int i = 0;i < SIZE;i++){
                for(int j = 0;j < SIZE;++j){
                    printf("%2.4f    ", P_Result[i][j]);
                }
                printf("\n");
            }
    } else {
         printf("对不起矩阵P_Result的结果为空！");
    }
    printf("控制率矩阵的结果：\n");
    for(int i = 0;i < SIZE;i++){
        printf("%2.4f    ", K[i]);
        printf("\n");
    }


}



