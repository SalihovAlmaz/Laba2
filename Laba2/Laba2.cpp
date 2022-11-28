
# define M_PI           3.14
#include <cmath>
#include <iostream>
#include <vector>
using namespace std;
double funcS(double x) {
    double t = sin(x);
    return t;
}
double fDerivate2(double x) {
    double t = -sin(x);
    return t;
}
/*
 Кубический сплайн
x  - текущий аргумент 
i - номер точки
М - массив вторых производных кубического сплайна, М[i] =  S3 ''[i]
X - массив аргументов 

*/
double S3(double x, int i, double h, vector<double> M, vector<double>X) {
    return M[i - 1] * (pow((X[i] - x), 3) - pow(h, 2) * (X[i] - x)) / (6 * h)
        + M[i] * (pow(x - X[i - 1], 3) - pow(h, 2) * (x - X[i - 1])) / (6 * h)
        + funcS(X[i - 1]) * (X[i] - x) / h
        + funcS(X[i]) * (x - X[i - 1]) / h;
}

/*
Метод прогонки
краевые условия S3 '' (a) = f''(a), S3''(b) = f''(b)
f - заданная функция
der2 - 2 производная заданной функции
f0 - начало отрезка
fn - конец отрезка
h -длина шага
n - кол-во частичных отрезков
X- массив аргументов
М - массив вторых производных кубическог сплайна

*/
vector<double> TMA(double(*funcS)(double),double(*der2)(double), double f0,double fn, double h, int n, vector<double> X) {
    // согласнов краевым условиям
    double ai = 2 * h / 3;
    double bi = h / 6;
    double ci = h / 6;
    double a0 = 1;
    double b0 = 0;
    double an = 1;
    double cn = 0;
    // правая часть СЛУ
     vector<double> d(n+1);
     d[0] = der2(f0); // f''(0) = -sin(0)
     d[n] = der2(fn); // f'' (Pi) = -sin(Pi)

     //Yi = (Yi+1 - Yi)/h - (Yi-Yi-1)/h
     for (int i = 1; i < n; i++) {
         d[i] = (funcS(X[i + 1]) - funcS(X[i])) / h - (funcS(X[i]) - funcS(X[i - 1])) / h;
     }

     //Прогоночные коэффициенты
     vector<double> l(n + 1);
     vector<double> u(n + 1);

     l[0] = -b0 / a0; //0
     u[0] = d[0] / a0; //0

     //Заполняем согласно рекуррентным формулам
     for (int i = 1; i < n; i++) {
         l[i] = -bi / (ai + ci * l[i - 1]);
         u[i] = (d[i] - ci * u[i - 1]) / (ai + ci * l[i - 1]);
     }

     // и крайнюю точку
     l[n] = -bi / (an + cn * l[n - 1]);
     u[n] = (d[n] - cn * u[n - 1]) / (an + cn * l[n - 1]);

     vector<double>M(n + 1);
     M[n] = u[n]; //т.к. bn = 0

     //Mi = li * Mi+1 + ui
     for (int i = n - 1; i >= 0; i--) {
         M[i] = l[i] * M[i + 1] + u[i];
     }
     return M;
}
int main()
{

    const double a = 0;
    const double b = M_PI;
    std::cout.precision(16);
    std::cout.setf(std::ios::fixed, std::ios::floatfield);
    vector<double> result;

    double prevMax = 0;
    double h;
    for (size_t n = 5; n < 5120; n *= 2)
    {
        //длина шага
        h = (b - a) / n;
        //Разбиваем на частичные отрезки
        vector<double> X{a};                    //111
        for (size_t i = 1; i < n; i++)
        {
            X.push_back(a + i * h);
        }
        X.push_back(b);

        //получаем вторые проивзодные кубического сплайна методо прогонки
        vector<double> M = TMA(funcS, fDerivate2, a, b, h, n, X);

        double deltaMax = 0;
        double ocenka = 0;
        for (size_t i = 1; i < n + 1; i++)
        {
            // |S3(x) - f(x)| в серединах частичных отрезков
            double s3 = S3(X[i - 1] + h / 2, i, h, M, X);
            ocenka = abs(s3 - funcS(X[i - 1] + h / 2));
            deltaMax = max(deltaMax, ocenka);                  
        }
        result.push_back(n);
        result.push_back(deltaMax);
        result.push_back(prevMax / 16);
        result.push_back(prevMax/deltaMax);
       
       
        
        prevMax = deltaMax;
    }
    
        cout << 'n' << "\t\t\t|\t" << "OcMax" << "\t\t\t|\t" << "Oc" << "\t\t\t|\t" << "K" << endl;
    
    for (size_t i = 0; i < result.size();)
    {
        for (size_t j = 0; j < 4; j++,i++) {
            cout << result[i] << "\t|\t";
        }
        cout << endl;
    }

}
