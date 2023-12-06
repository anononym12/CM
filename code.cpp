// BISECTION METHOD

#include <iostream>
#include <cmath>
using namespace std;
#define tol 0.000001
double fx(double x) {
    return x * x * x - 8 * x - 3;
}
double bisection(double a, double b) {
    double mid, approx;
    int iteration = 0;
    while (fabs(b - a) > tol) {
        mid = (a + b) / 2.0;
        approx = fx(mid);
        cout << "Iteration " << ++iteration << ": Approximation = " << mid << endl;

        if (approx == 0) {
            break;
        } else if (approx * fx(a) < 0) {
            b = mid;
        } else {
            a = mid;
        }
    }
    return approx == 0 ? mid : (a + b) / 2.0;
}
int main() {
    double a, b;
    cout << "Enter initial interval [a, b]: ";
    cin >> a >> b;
    double root = bisection(a, b);
    cout << "\nApproximated Root: " << root << endl;
    return 0;
}


// NEWTON RAPHSON METHOD
#include <iostream>
#include <cmath>

using namespace std;

#define tol 0.000001
#define maxIterations 1000

double fx(double x) {
    return x * x * x - 8 * x - 3;
}

double derivative(double x) {
    return 3 * x * x - 8;
}

double newtonRaphson(double initialGuess) {
    double x = initialGuess;
    double approx;

    int iteration = 0;
    do {
        approx = fx(x) / derivative(x);
        x = x - approx;

        cout << "Iteration " << ++iteration << ": Approximation = " << x << endl;

        if (fabs(approx) <= tol) {
            break; // Convergence reached
        }
    } while (iteration < maxIterations);

    return x;
}

int main() {
    double initialGuess;

    cout << "Enter initial guess: ";
    cin >> initialGuess;

    double root = newtonRaphson(initialGuess);

    cout << "\nApproximated Root: " << root << endl;

    return 0;
}

// SECANT METHOD
#include <iostream>
#include <cmath>

using namespace std;

#define tol 0.00001
#define maxIterations 1000

double f(double x) {
    return pow(x, 3) + 10 * x - 20;
}

double secantMethod(double x0, double x1) {
    double xn;
    for (int n = 0; n < maxIterations; ++n) {
        xn = x1 - (f(x1) * (x1 - x0)) / (f(x1) - f(x0));

        cout << "Iteration (" << n << "): x"<<n << "= " << xn << endl;

        if (fabs(xn - x1) < tol) {
            cout << "Convergence achieved." << endl;
            return xn;
        }

        x0 = x1;
        x1 = xn;
    }

    cout << "Maximum number of iterations reached. Convergence not achieved." << endl;
    return xn;
}

int main() {
    double x0, x1;

    cout << "Enter initial guess x0: ";
    cin >> x0;

    cout << "Enter initial guess x1: ";
    cin >> x1;

    double root = secantMethod(x0, x1);

    cout << "\nApproximated Root: " << root << endl;

    return 0;
}

// POWER METHOD

#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

const int maxIterations = 1000;
const double tol = 0.000001;

vector<double> matrixVectorMultiply(const vector<vector<double>>& A, const vector<double>& v) {
    int n = A.size();
    vector<double> result(n, 0.0);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            result[i] += A[i][j] * v[j];
        }
    }

    return result;
}

void powerMethod(const vector<vector<double>>& A, double& largestEigenvalue, double& smallestEigenvalue) {
    int n = A.size();
    vector<double> v(n, 1.0);

    for (int iteration = 0; iteration < maxIterations; ++iteration) {
        vector<double> Av = matrixVectorMultiply(A, v);
        double maxElement = 0.0;

        for (int i = 0; i < n; ++i) {
            maxElement = max(maxElement, abs(Av[i]));
        }

        for (int i = 0; i < n; ++i) {
            v[i] = Av[i] / maxElement;
        }

        largestEigenvalue = maxElement;
        smallestEigenvalue = 1.0 / maxElement;

        if (abs(Av[0] - largestEigenvalue * v[0]) < tol) {
            break;
        }
    }
}

int main() {
    int n;

    cout << "Enter the size of the square matrix: ";
    cin >> n;

    vector<vector<double>> matrix(n, vector<double>(n, 0.0));

    cout << "Enter the elements of the matrix:" << endl;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cin >> matrix[i][j];
        }
    }

    double largestEigenvalue, smallestEigenvalue;

    powerMethod(matrix, largestEigenvalue, smallestEigenvalue);

    cout << "Largest Eigenvalue: " << largestEigenvalue << endl;
    cout << "Smallest Eigenvalue: " << smallestEigenvalue << endl;

    return 0;
}


// RUNGE-METHOD

#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

double equation(double t, double x) {
    return 100 * (sin(t) - x);
}

void rungeKutta(double t0, double x0, double tn, double h) {
    int numSteps = static_cast<int>((tn - t0) / h) + 1;
    vector<double> t(numSteps), x(numSteps);

    t[0] = t0;
    x[0] = x0;

    for (int i = 1; i < numSteps; ++i) {
        double k1 = h * equation(t[i - 1], x[i - 1]);
        double k2 = h * equation(t[i - 1] + h / 2, x[i - 1] + k1 / 2);
        double k3 = h * equation(t[i - 1] + h / 2, x[i - 1] + k2 / 2);
        double k4 = h * equation(t[i - 1] + h, x[i - 1] + k3);

        t[i] = t[i - 1] + h;
        x[i] = x[i - 1] + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    }

    cout << "Step Size: " << h << "\tLast Iteration: " << numSteps - 1 << "\tResult: " << x.back() << endl;
}

int main() {
    double t0 = 0.0;
    double x0 = 0.0;
    double tn = 3.0;

    rungeKutta(t0, x0, tn, 0.015);
    rungeKutta(t0, x0, tn, 0.020);
    rungeKutta(t0, x0, tn, 0.025);
    rungeKutta(t0, x0, tn, 0.030);

    return 0;
}


// LAGRANGE INTERPOLATION
#include <iostream>
#include <vector>

using namespace std;

double lagrangeInterpolation(const vector<int>& years, const vector<double>& rainfall, int n, int x) {
    double result = 0.0;

    for (int i = 0; i < n; i++) {
        double term = rainfall[i];
        for (int j = 0; j < n; j++) {
            if (j != i) {
                term = term * (x - years[j]) / (years[i] - years[j]);
            }
        }
        result += term;
    }

    return result;
}

int main() {
    int numYears;
    cout << "Enter the number of years of rainfall data: ";
    cin >> numYears;

    vector<int> years(numYears);
    vector<double> rainfall(numYears);

    cout << "Enter the annual rainfall data:" << endl;
    for (int i = 0; i < numYears; i++) {
        cout << "Year " << i + 1 << ": ";
        cin >> years[i];
        cout << "Rainfall for Year " << i + 1 << ": ";
        cin >> rainfall[i];
    }

    int currentYear;
    cout << "Enter the current year for rainfall prediction: ";
    cin >> currentYear;

    double predictedRainfall = lagrangeInterpolation(years, rainfall, numYears, currentYear);

    cout << "Predicted rainfall for " << currentYear << ": " << predictedRainfall << " mm" << endl;

    return 0;
}


// TRAPEZOIDAL RULE
#include<iostream>
#include<cmath>
#include<iomanip>
using namespace std;

double error = 0.00005;

double f(double x) {
    return exp(-x * x);
}

double I(double a, double b, double h, int n) {
    double val = 0;
    for (int i = 1; i < n; i++) {
        val += f(a + i * h);
    }
    return h / 2 * (f(a) + f(b) + 2 * (val));
}

int main() {
    cout << "Using COMPOSITE TRAPEZOIDAL RULE:" << endl;
    int n = 1; // Start with a non-zero number of intervals
    double a = 0, b = 1, Integral = 0;
    double h = (b - a) / n;
    cout << "interval \t h \t\t Integrated value \t error" << endl;
    
    do {
        Integral = I(a, b, h, n);
        cout << fixed << setprecision(5) << n << " \t\t " << h << " \t " << Integral << " \t\t " << fabs(Integral - 0.746824132812) << endl;
        n++;
        h = (b - a) / n;
    } while (fabs(Integral - 0.746824132812) > error);

    cout << "\nNumber of intervals required to evaluate is: " << n - 1 << endl << endl;
    return 0;
}



// FIBONACCI 
#include <iostream>
#include <cmath>

using namespace std;

double f(double x) {
    return (x * x + sin(53 * x));
}

int main() {
    cout << "\nMinima of Function using FIBONACCI SEARCH ALGORITHM METHOD : \n";
    double a = -0.2, b = 0.2;
    double x, y;
    double D, d;
    int N = ceil(log((b - a) / (0.0001 * 1.618)) / log(1.618)) + 1;
    cout << "\nNumber of steps (N) = " << N - 1 << endl;

    int fib[N];
    fib[0] = 0;
    fib[1] = 1;
    for (int i = 2; i < N; i++) {
        fib[i] = fib[i - 1] + fib[i - 2];
    }
    cout << "\nFibonacci Numbers used are : " << endl;
    for (int i = 0; i < N; i++) {
        cout << fib[i] << " ";
    }
    cout << endl
         << endl;

    cout << "k \t a \t\t b \t\t x \t\t y\n";
    int k;
    for (k = 2; k < N; k++) {
        D = (1.0 * fib[k - 2] / fib[k]) * (b - a);
        x = a + D;
        y = b - D;
        cout << fixed << x << " \t " << a << " \t " << b << " \t " << x << " \t " << y << endl;
        if (f(x) >= f(y))
            a = x;
        else
            b = y;
    }
    d = (y - x) / 4;
    x = (a + b) / 2 - 2 * d;
    y = (a + b) / 2 + 2 * d;
    cout << fixed << x << " \t " << a << " \t " << b << " \t " << x << " \t " << y << endl;
    if (f(x) >= f(y))
        a = x;
    else
        b = y;

    double xmin = (a + b) / 2;
    double minima = f(xmin);

    cout << "Minima is at (a+b)/2 = " << xmin << endl;
    cout << "\nMinima of the function near the origin at " << xmin << " is " << minima;

    cout << endl
         << endl;
    return 0;
}



