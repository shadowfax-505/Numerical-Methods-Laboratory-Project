# Numerical Methods
### Introduction
The project is a classic implementation of the 19 numerical methods we learned in our numerical methods labratory. Each methods are use to solve a particular type of equation . Each methods have their own respective features . We built the project so that the methods take input from a .txt file and gives output to another .txt file. All methods are well structured under their respective types.Such as:

### 1. Solution of Non-Linear Equations (4 methods)
* Bisection Method
* False Position Method
* Newton-Raphson Method
* Secant Method
### 2. Solution of Linear Equations (4 methods)
* Gauss Elimination
* Gauss-Jordan Elimination
* LU Decomposition (Doolittle)
* Matrix Inversion
### 3. Differential Equation Solving (1 method)
* Runge-Kutta 4th Order Method
### 4. Interpolation Methods (3 methods)
* Newton Forward Interpolation
* Newton Backward Interpolation
* Newton Divided Difference Interpolation
### 5. Numerical Differentiation (2 methods)
* Differentiation by Forward Interpolation
* Differentiation by Backward Interpolation
### 6. Curve Fitting & Regression (3 methods)
* Linear Regression
* Polynomial Regression
* Transcendental Regression
### 7. Numerical Integration (2 methods)
* Simpson's 1/3 Rule
* Simpson's 3/8 Rule



### Project Structure

```
README.md/
├── Project Introduction
├── README.md Structure
└── Table of Contents
├── Solution of Non-Linear Equations
    ├── Bisection
    ├── False Position
    ├── Secant
    ├── Newton Raphson

├── Solution of Linear Equations
    ├── Gauss Elimination
    ├── Gauss Jordan Elimination
    ├── LU Decomposition
    ├── Matrix Inversion

├── Solution of Differential Equations
    ├── Runge Kutta 4th Order

├── Numerical Integration
    ├── Simpson's 1/3 Rule
    └── Simpson's 3/8 Rule

├── Interpolation Methods
    ├── Newton Forward Interpolation
    ├── Newton Backward Interpolation
    ├── Newton Divided Difference Interpolation

├── Numerical Differentiation
    ├── Differentiation by Forward Interpolation
    ├── Differentiation by Backward Interpolation

├── Curve Fitting: Regression
    ├── Linear Regression
    ├── Polynomial Regression
    ├── Transcendental Regression


```


### Table of Contents

- [Solution of Non-Linear Equations](#solution-of-non-linear-equations)
  - [Bisection Method](#bisection-method)
    - [Theory](#bisection-theory)
    - [Code](#bisection-code)
    - [Input](#bisection-input)
    - [Output](#bisection-output)
  - [False Position Method](#false-position-method)
    - [Theory](#false-position-theory)
    - [Code](#false-position-code)
    - [Input](#false-position-input)
    - [Output](#false-position-output)
  - [Secant Method](#secant-method)
    - [Theory](#secant-theory)
    - [Code](#secant-code)
    - [Input](#secant-input)
    - [Output](#secant-output)
  - [Newton Raphson Method](#newton-raphson-method)
    - [Theory](#newton-raphson-theory)
    - [Code](#newton-raphson-code)
    - [Input](#newton-raphson-input)
    - [Output](#newton-raphson-output)

- [Solution of Linear Equations](#solution-of-linear-equations)
  - [Gauss Elimination Method](#gauss-elimination-method)
    - [Theory](#gauss-elimination-theory)
    - [Code](#gauss-elimination-code)
    - [Input](#gauss-elimination-input)
    - [Output](#gauss-elimination-output)
  - [Gauss Jordan Elimination Method](#gauss-jordan-elimination-method)
    - [Theory](#gauss-jordan-theory)
    - [Code](#gauss-jordan-code)
    - [Input](#gauss-jordan-input)
    - [Output](#gauss-jordan-output)
  - [LU Decomposition Method](#lu-decomposition-method)
    - [Theory](#lu-decomposition-theory)
    - [Code](#lu-decomposition-code)
    - [Input](#lu-decomposition-input)
    - [Output](#lu-decomposition-output)
  - [Matrix Inversion](#matrix-inversion-method)
    - [Theory](#matrix-inversion-theory)
    - [Code](#matrix-inversion-code)
    - [Input](#matrix-inversion-input)
    - [Output](#matrix-inversion-output)

- [Solution of Differential Equations](#solution-of-differential-equations)
  - [Runge-Kutta 4th Order Method](#runge-kutta-4th-order-method)
    - [Theory](#runge-kutta-4th-order-theory)
    - [Code](#runge-kutta-4th-order-code)
    - [Input](#runge-kutta-4th-order-input)
    - [Output](#runge-kutta-4th-order-output)

- [Numerical Integration](#numerical-integration)
  - [Simpson's 1/3 Rule](#simpsons-13-rule)
    - [Theory](#simpsons-13-theory)
    - [Code](#simpsons-13-code)
    - [Input](#simpsons-13-input)
    - [Output](#simpsons-13-output)
  - [Simpson's 3/8 Rule](#simpsons-38-rule)
    - [Theory](#simpsons-38-theory)
    - [Code](#simpsons-38-code)
    - [Input](#simpsons-38-input)
    - [Output](#simpsons-38-output)

- [Interpolation Methods](#interpolation-methods)
  - [Newton Forward Interpolation](#newton-forward-interpolation-method)
    - [Theory](#newton-forward-theory)
    - [Code](#newton-forward-code)
    - [Input](#newton-forward-input)
    - [Output](#newton-forward-output)
  - [Newton Backward Interpolation](#newton-backward-interpolation-method)
    - [Theory](#newton-backward-theory)
    - [Code](#newton-backward-code)
    - [Input](#newton-backward-input)
    - [Output](#newton-backward-output)
  - [Newton Divided Difference Interpolation](#newton-divided-difference-interpolation-method)
    - [Theory](#newton-divided-difference-theory)
    - [Code](#newton-divided-difference-code)
    - [Input](#newton-divided-difference-input)
    - [Output](#newton-divided-difference-output)

- [Numerical Differentiation](#numerical-differentiation)
  - [Differentiation by Forward Interpolation](#differentiation-forward-interpolation-method)
    - [Theory](#differentiation-forward-theory)
    - [Code](#differentiation-forward-code)
    - [Input](#differentiation-forward-input)
    - [Output](#differentiation-forward-output)
  - [Differentiation by Backward Interpolation](#differentiation-backward-interpolation-method)
    - [Theory](#differentiation-backward-theory)
    - [Code](#differentiation-backward-code)
    - [Input](#differentiation-backward-input)
    - [Output](#differentiation-backward-output)

- [Curve Fitting: Regression](#curve-fitting-regression)
  - [Linear Regression](#linear-regression-method)
    - [Theory](#linear-regression-theory)
    - [Code](#linear-regression-code)
    - [Input](#linear-regression-input)
    - [Output](#linear-regression-output)
  - [Polynomial Regression](#polynomial-regression-method)
    - [Theory](#polynomial-regression-theory)
    - [Code](#polynomial-regression-code)
    - [Input](#polynomial-regression-input)
    - [Output](#polynomial-regression-output)
  - [Transcendental Regression](#transcendental-regression-method)
    - [Theory](#transcendental-regression-theory)
    - [Code](#transcendental-regression-code)
    - [Input](#transcendental-regression-input)
    - [Output](#transcendental-regression-output)


    
   
    
---


## Solution of Non-Linear Equations

### Bisection Method

#### Bisection Theory

The Bisection Method is a numerical method which is used to find the real roots of a polynomial function. It works by dividing an interval into two halves and selecting the sub-interval in which the root lies. The method is based on the Intermediate Value Theorem, which guarantees the existence of a root within an interval where the function values at the endpoints have opposite signs. This method gradually reduces the interval containing the root by taking the average (midpoint) of the interval endpoints. Although simple and reliable, the bisection method converges relatively slowly compared to other numerical methods.

Mathematical Principle
Let f(x) be a continuous function defined on a closed interval [a, b], such that:
f(a) * f(b) < 0
Then, according to the Intermediate Value Theorem, there exists at least one value x in (a, b) for which: f(x) = 0

Bisection Method Algorithm
To find the root of a continuous function f(x), follow these steps:

1. Choose two initial points a and b such that a < b and f(a) * f(b) < 0
2. Compute the midpoint:
   t = (a + b) / 2
3. If f(t) = 0, then t is the root
4. Otherwise:
   - If f(t) * f(a) < 0, the root lies in the interval [a, t]
   - Else if f(t) * f(b) < 0, the root lies in the interval [t, b]
5. Repeat the above steps until the interval becomes sufficiently small or the desired accuracy is achieved

The bisection method is an approximation technique that continues dividing the interval until the root is obtained within a set error margin.

Features
- Supports multiple test cases
- Dynamic evaluation of polynomial functions
- Automatic generation of function strings
- Eliminates duplicate roots
- Well-formatted output with precision control
- Default error tolerance: 0.001

#### Bisection Code
```cpp
#include <bits/stdc++.h>
using namespace std;

vector<double> coef;
int deg;

double fx(double x)
{
    double s = 0;
    for (int i = 0; i <= deg; i++)
        s += coef[i] * pow(x, deg - i);
    return s;
}

string funcStr()
{
    stringstream ss;
    ss << "f(x)=";
    bool first = true;

    for (int i = 0; i <= deg; i++)
    {
        double c = coef[i];
        if (fabs(c) < 1e-10) continue;

        if (!first)
            ss << (c > 0 ? "+" : "-");
        else if (c < 0)
            ss << "-";

        double v = fabs(c);
        int p = deg - i;

        if (p == 0)
            ss << v;
        else if (p == 1)
            ss << (fabs(v - 1.0) < 1e-10 ? "x" : to_string(v) + "x");
        else
            ss << (fabs(v - 1.0) < 1e-10 ? "x^" : to_string(v) + "x^") << p;

        first = false;
    }
    return ss.str();
}

double limitCalc()
{
    if (deg < 2) return 10.0;
    double a = coef[0], b = coef[1], c = coef[2];
    if (fabs(a) < 1e-10) return 10.0;
    return sqrt((b / a) * (b / a) - 2 * (c / a));
}

void uniq(vector<double> &r, double e)
{
    sort(r.begin(), r.end());
    r.erase(unique(r.begin(), r.end(),
        [&](double x, double y)
        {
            return fabs(x - y) < e;
        }), r.end());
}

void solve(int tc, ofstream &out)
{
    vector<double> roots;
    double eps = 0.001;
    double lim = limitCalc();

    for (double i = -lim; i <= lim; i += 0.5)
    {
        double a = i, b = i + 0.5;

        if (fx(a) * fx(b) < 0)
        {
            while (fabs(b - a) > eps)
            {
                double m = (a + b) / 2;
                if (fx(a) * fx(m) < 0)
                    b = m;
                else
                    a = m;
            }
            roots.push_back((a + b) / 2);
        }
    }

    uniq(roots, eps);

    out << "TestCase " << tc << "\n";
    out << "Function:" << funcStr() << "\n";
    out << "Degree:" << deg << "\n";
    out << "Tolerance:" << eps << "\n";
    out << "Roots:\n";

    if (roots.empty())
        out << "None\n";
    else
        for (double x : roots)
            out << fixed << setprecision(6) << x << "\n";

    out << "\n";
}

int main()
{
    ifstream in("Input.txt");
    ofstream out("Output.txt");

    int T;
    in >> T;
    out << "BisectionMethod\n";
    out << "Cases:" << T << "\n\n";

    for (int t = 1; t <= T; t++)
    {
        in >> deg;
        coef.assign(deg + 1, 0);

        for (int i = 0; i <= deg; i++)
            in >> coef[i];

        solve(t, out);
    }

    return 0;
}

```

#### Bisection Input
```
3
3
1 -6 11 -6
2
1 -4 4
5
1 -3 -7 27 -10 -8

```

#### Bisection Output
```
BisectionMethod
Cases:3

TestCase:1
Function:f(x)=x^3-6x^2+11x-6
Degree:3
Tolerance:0.001
Roots:
1.000000
2.000000
3.000000

TestCase:2
Function:f(x)=x^2-4x+4
Degree:2
Tolerance:0.001
Roots:
2.000000

TestCase:3
Function:f(x)=x^5-3x^4-7x^3+27x^2-10x-8
Degree:5
Tolerance:0.001
Roots:
-2.000244
0.999756
2.000244

```
[Back to Table of Contents](#table-of-contents)

---

### False Position Method

#### False Position Theory

The False Position Method, also known as the Regula Falsi Method, is a numerical technique used to find the roots of a continuous function. Like the Bisection Method, it requires an interval [a, b] where the function changes sign, i.e., f(a) * f(b) < 0.
The False Position Method uses a straight line connecting the points (a, f(a)) and (b, f(b)) and finds the point where this line intersects the x-axis. This intersection is used as a better approximation of the root.
The method iteratively updates the interval until the root is found with the desired accuracy. It often converges faster than the Bisection Method, especially when the function is nearly linear near the root.

Mathematical Principle
Given a continuous function f(x) on [a, b] such that:
f(a) * f(b) < 0
The point c (approximation of the root) is calculated as:
c = b - (f(b) * (a - b)) / (f(a) - f(b))
Then, depending on the sign of f(c):
- If f(a) * f(c) < 0, the root lies in [a, c]
- Else if f(b) * f(c) < 0, the root lies in [c, b]
Repeat the process until |f(c)| is smaller than the desired tolerance.

Algorithm (False Position Method)
1. Choose two initial points a and b such that f(a) * f(b) < 0
2. Compute the intersection point of the line joining (a, f(a)) and (b, f(b)):
   c = b - (f(b) * (a - b)) / (f(a) - f(b))
3. If f(c) = 0, then c is the root.
4. Otherwise, replace either a or b depending on the sign of f(c) to form a new interval where a root exists.
5. Repeat steps 2–4 until the root is approximated within the desired error tolerance.

Features
- Supports multiple test cases
- Faster convergence than Bisection in many cases
- Dynamic function evaluation
- Precision control in output
- Default error tolerance: 0.001


#### False Position Code
```cpp
#include <bits/stdc++.h>
using namespace std;

double evalPoly(const vector<double>& coef, double x)
{
    double val = 0;
    int deg = coef.size() - 1;
    for(int i = 0; i <= deg; i++)
        val += coef[i] * pow(x, deg - i);
    return val;
}

string polyStr(const vector<double>& coef)
{
    stringstream ss;
    int deg = coef.size() - 1;
    bool first = true;

    for(int i = 0; i <= deg; i++)
    {
        double c = coef[i];
        int p = deg - i;
        if(fabs(c) < 1e-10) continue;

        if(!first) ss << (c > 0 ? " + " : " - ");
        else if(c < 0) ss << "-";

        double v = fabs(c);

        if(p == 0) ss << fixed << setprecision(4) << v;
        else if(p == 1) ss << fixed << setprecision(4) << v << "*x";
        else ss << fixed << setprecision(4) << v << "*x^" << p;

        first = false;
    }
    return ss.str();
}

double searchLimit(const vector<double>& coef)
{
    if(coef.size() < 3) return 10.0;
    double a = coef[0], b = coef[1], c = coef[2];
    return sqrt(abs((b / a) * (b / a) - 2 * (c / a)));
}

int main()
{
    ifstream fin("Input.txt");
    ofstream fout("Output.txt");

    if(!fin.is_open())
    {
        cerr << "Error: cannot open Input.txt\n";
        return 1;
    }

    int T;
    fin >> T;

    for(int t = 1; t <= T; t++)
    {
        int deg;
        fin >> deg;
        vector<double> coef(deg + 1);
        for(int i = 0; i <= deg; i++) fin >> coef[i];

        double tol;
        fin >> tol;

        cout << "TestCase " << t << "\n";
        fout << "TestCase " << t << "\n";

        string pstr = polyStr(coef);
        cout << "Polynomial: f(x)=" << pstr << "\n";
        fout << "Polynomial: f(x)=" << pstr << "\n";

        cout << "Tolerance: " << scientific << tol << "\n";
        fout << "Tolerance: " << scientific << tol << "\n";

        double lim = min(searchLimit(coef), 100.0);
        cout << "Search range: [" << -lim << ", " << lim << "]\n";
        fout << "Search range: [" << -lim << ", " << lim << "]\n";

        vector<double> roots;
        int rootCount = 0;

        for(double x = -lim; x < lim; x += 0.5)
        {
            double xl = x, xr = x + 0.5;
            double fl = evalPoly(coef, xl);
            double fr = evalPoly(coef, xr);

            if(fabs(fl) < tol) { roots.push_back(xl); continue; }
            if(fabs(fr) < tol) { roots.push_back(xr); continue; }

            if(fl * fr < 0)
            {
                int iter = 0;
                double xm;
                while(fabs(xr - xl) > tol)
                {
                    iter++;
                    xm = (xl * fr - xr * fl) / (fr - fl);
                    double fm = evalPoly(coef, xm);

                    if(fabs(fm) < tol) break;

                    if(fm * fl < 0) { xr = xm; fr = fm; }
                    else { xl = xm; fl = fm; }
                }

                roots.push_back(xm);
                rootCount++;
                cout << "Root " << rootCount << " in [" << fixed << setprecision(4) << xl << ", " << xr << "] = " << xm << " (iter: " << iter << ")\n";
                fout << "Root " << rootCount << " in [" << fixed << setprecision(4) << xl << ", " << xr << "] = " << xm << " (iter: " << iter << ")\n";
            }
        }

        sort(roots.begin(), roots.end());
        roots.erase(unique(roots.begin(), roots.end(), [tol](double a, double b){ return fabs(a-b) < tol; }), roots.end());

        cout << "All Real Roots:\n";
        fout << "All Real Roots:\n";

        if(roots.empty()) { cout << "None\n"; fout << "None\n"; }
        else
        {
            for(size_t i = 0; i < roots.size(); i++)
            {
                cout << "Root " << i+1 << ": " << fixed << setprecision(6) << roots[i] << "\n";
                fout << "Root " << i+1 << ": " << fixed << setprecision(6) << roots[i] << "\n";
            }
        }

        cout << "\n";
        fout << "\n";
    }

    fin.close();
    fout.close();
    return 0;
}

```

#### False Position Input
```
3
3
1 -6 11 -6
0.001
2
1 -4 4
0.001
5
1 -3 -7 27 -10 -8
0.001

```

#### False Position Output
```
TestCase 1
Polynomial: f(x)=1.0000*x^3-6.0000*x^2+11.0000*x-6.0000
Tolerance: 1.000000e-03
Search range: [-3.0, 3.0]
Root 1 in [0.5000, 1.0000] = 1.000000 (iter: 11)
Root 2 in [1.5000, 2.0000] = 2.000000 (iter: 12)
Root 3 in [2.5000, 3.0000] = 3.000000 (iter: 12)
All Real Roots:
Root 1: 1.000000
Root 2: 2.000000
Root 3: 3.000000

TestCase 2
Polynomial: f(x)=1.0000*x^2-4.0000*x+4.0000
Tolerance: 1.000000e-03
Search range: [-10.0, 10.0]
Root 1 in [1.5000, 2.0000] = 2.000000 (iter: 13)
All Real Roots:
Root 1: 2.000000

TestCase 3
Polynomial: f(x)=1.0000*x^5-3.0000*x^4-7.0000*x^3+27.0000*x^2-10.0000*x-8.0000
Tolerance: 1.000000e-03
Search range: [-5.47723, 5.47723]
Root 1 in [-2.0000, -1.5000] = -2.000244 (iter: 14)
Root 2 in [0.5000, 1.0000] = 0.999756 (iter: 15)
Root 3 in [1.5000, 2.0000] = 2.000244 (iter: 15)
All Real Roots:
Root 1: -2.000244
Root 2: 0.999756
Root 3: 2.000244

```
[Back to Table of Contents](#table-of-contents)

---

### Secant Method

#### Secant Theory

The Secant Method is a numerical method used to find the roots of a continuous function. It is similar to the Newton-Raphson method but does not require the computation of derivatives. Instead, it uses a line through two points on the function curve to approximate the root.
The method iteratively updates the approximation of the root using the two most recent estimates. The Secant Method is generally faster than the Bisection Method but may be less stable.

Mathematical Principle
Given a continuous function f(x) and two initial approximations x0 and x1, the next approximation x2 is computed using the formula:
x2 = x1 - f(x1) * ((x1 - x0) / (f(x1) - f(x0)))

Then, the process is repeated using the latest two approximations:
x_{n+1} = x_n - f(x_n) * ((x_n - x_{n-1}) / (f(x_n) - f(x_{n-1}))) 

This continues until the value of f(x_{n+1}) is sufficiently close to zero or the difference between successive approximations is smaller than the desired tolerance.

Algorithm (Secant Method)
1. Choose two initial points x0 and x1 near the root.
2. Compute the next approximation using:
   x2 = x1 - f(x1) * ((x1 - x0) / (f(x1) - f(x0)))
3. Check if |f(x2)| < tolerance or |x2 - x1| < tolerance.
4. If not, set x0 = x1 and x1 = x2, then repeat step 2.
5. Continue until the root is approximated within the desired error tolerance.

Features
- Does not require derivative calculation
- Faster convergence than Bisection in most cases
- Requires good initial guesses
- Can fail to converge if function is not well-behaved
- Supports multiple test cases
- Default error tolerance: 0.001


#### Secant Code
```cpp
#include <bits/stdc++.h>
using namespace std;

vector<double> foundRoots;
ofstream fout;

double evalPoly(double x, const vector<double>& coef)
{
    double sum = 0;
    for(int i = 0; i < coef.size(); i++)
        sum += coef[i] * pow(x, i);
    return sum;
}

void printPolynomial(const vector<double>& coef)
{
    bool firstTerm = true;
    int n = coef.size();

    for(int i = n-1; i >= 0; i--)
    {
        if(fabs(coef[i]) < 1e-12) continue;

        if(!firstTerm)
        {
            if(coef[i] < 0) cout << " - ", fout << " - ";
            else cout << " + ", fout << " + ";
        }
        else
        {
            if(coef[i] < 0) cout << "-", fout << "-";
            firstTerm = false;
        }

        double c = fabs(coef[i]);

        if(i == 0) cout << c, fout << c;
        else if(i == 1) cout << c << "x", fout << c << "x";
        else cout << c << "x^" << i, fout << "x^" << i;
    }

    cout << endl;
    fout << endl;
}

void secantMethod(double x0, double x1, double tol, const vector<double>& coef)
{
    double xNext;
    int iterations = 0;

    while(fabs(x1 - x0) >= tol && fabs(evalPoly(x1, coef)) >= tol)
    {
        double f0 = evalPoly(x0, coef);
        double f1 = evalPoly(x1, coef);

        if(fabs(f1 - f0) < 1e-12) return;

        xNext = (x0 * f1 - x1 * f0) / (f1 - f0);
        iterations++;

        x0 = x1;
        x1 = xNext;
    }

    if(fabs(evalPoly(xNext, coef)) > tol) return;

    for(double r : foundRoots)
        if(fabs(r - xNext) < tol) return;

    foundRoots.push_back(xNext);

    cout << fixed << setprecision(6);
    fout << fixed << setprecision(6);

    cout << "Interval: [" << x0 << ", " << x1 << "] Root: " << xNext
         << " Iterations: " << iterations << endl;
    fout << "Interval: [" << x0 << ", " << x1 << "] Root: " << xNext
         << " Iterations: " << iterations << endl;
}

int main()
{
    ifstream fin("Input.txt");
    fout.open("Output.txt");

    if(!fin.is_open())
    {
        cout << "Cannot open Input.txt" << endl;
        return 0;
    }

    int testCases;
    fin >> testCases;

    while(testCases--)
    {
        foundRoots.clear();

        int deg;
        fin >> deg;
        vector<double> coef(deg + 1);
        for(double &x : coef) fin >> x;
        reverse(coef.begin(), coef.end());

        cout << "\nPolynomial: ";
        fout << "\nPolynomial: ";
        printPolynomial(coef);

        double xmax = 0;
        for(int i = 0; i <= deg; i++)
            xmax = max(xmax, 1.0 + fabs(coef[i]/coef[deg]));

        cout << "Root bound: " << xmax << endl;
        fout << "Root bound: " << xmax << endl;

        double step = 0.45;
        double tol = 1e-3;

        for(double i = -xmax; i <= xmax; i += step)
            secantMethod(i, i + step, tol, coef);
    }

    fin.close();
    fout.close();
    return 0;
}

```

#### Secant Input
```
2
3
1 -6 11 -6
5
1 -3 -7 27 -10 -8

```

#### Secant Output
```
Polynomial: 1x^3 - 6x^2 + 11x - 6
Root bound: 3
Interval: [0.5000, 1.0000] Root: 1.000000 Iterations: 6
Interval: [1.5000, 2.0000] Root: 2.000000 Iterations: 7
Interval: [2.5000, 3.0000] Root: 3.000000 Iterations: 7

Polynomial: 1x^5 - 3x^4 - 7x^3 + 27x^2 - 10x - 8
Root bound: 5.47723
Interval: [-2.0000, -1.5000] Root: -2.000244 Iterations: 8
Interval: [0.5000, 1.0000] Root: 0.999756 Iterations: 9
Interval: [1.5000, 2.0000] Root: 2.000244 Iterations: 9

```
[Back to Table of Contents](#table-of-contents)

---

### Newton Raphson Method

#### Newton Raphson Theory

The Newton-Raphson Method is a widely used numerical technique to find the roots of a continuous and differentiable function. It uses the derivative of the function to iteratively approximate the root. This method generally converges faster than the Bisection and Secant methods when the initial guess is sufficiently close to the actual root.

Mathematical Principle
Given a function f(x) and its derivative f'(x), if x_n is the current approximation of the root, the next approximation x_{n+1} is calculated as:
x_{n+1} = x_n - f(x_n) / f'(x_n)

This formula is applied repeatedly until the value of f(x_{n+1}) is close enough to zero or the difference between successive approximations is smaller than a specified tolerance.

Algorithm (Newton-Raphson Method)
1. Choose an initial guess x0 near the root.
2. Compute the next approximation using:
   x_{n+1} = x_n - f(x_n) / f'(x_n)
3. Check if |f(x_{n+1})| < tolerance or |x_{n+1} - x_n| < tolerance.
4. If not, set x_n = x_{n+1} and repeat step 2.
5. Continue until the root is approximated within the desired error tolerance.

Features
- Fast convergence if the initial guess is close to the root
- Requires derivative of the function
- May fail to converge if the initial guess is poor or function is not well-behaved
- Supports multiple test cases
- Default error tolerance: 0.001


#### Newton Raphson Code
```cpp
#include <bits/stdc++.h>
using namespace std;

double polyVal(const vector<double>& coef, double x)
{
    double sum = 0;
    int n = coef.size();
    for(int i = 0; i < n; i++)
        sum += coef[i] * pow(x, n-1-i);
    return sum;
}

double polyDeriv(const vector<double>& coef, double x)
{
    double sum = 0;
    int n = coef.size();
    for(int i = 0; i < n-1; i++)
        sum += (n-1-i) * coef[i] * pow(x, n-2-i);
    return sum;
}


string polyToStr(const vector<double>& coef)
{
    stringstream ss;
    int n = coef.size();
    bool first = true;

    for(int i = 0; i < n; i++)
    {
        double v = coef[i];
        int p = n-1-i;
        if(fabs(v) < 1e-12) continue;

        if(!first && v >= 0) ss << " + ";
        else if(v < 0) { ss << " - "; v = -v; }

        if(p == 0) ss << v;
        else if(p == 1) ss << v << "x";
        else ss << v << "x^" << p;

        first = false;
    }
    return ss.str();
}

int main()
{
    ifstream fin("Input.txt");
    ofstream fout("Output.txt");

    if(!fin.is_open())
    {
        cerr << "Error: cannot open Input.txt\n";
        return 1;
    }

    int T;
    fin >> T;

    for(int tc = 1; tc <= T; tc++)
    {
        int deg;
        fin >> deg;
        vector<double> coef(deg+1);
        for(int i = 0; i <= deg; i++) fin >> coef[i];

        double start, end, tol;
        fin >> start >> end >> tol;

        int maxIter;
        fin >> maxIter;

        cout << "TestCase " << tc << "\n";
        fout << "TestCase " << tc << "\n";

        string pstr = polyToStr(coef);
        cout << "Polynomial: f(x)=" << pstr << "\n";
        fout << "Polynomial: f(x)=" << pstr << "\n";

        cout << "Search interval: [" << start << ", " << end << "]\n";
        fout << "Search interval: [" << start << ", " << end << "]\n";

        cout << "Tolerance: " << scientific << tol << "\n";
        fout << "Tolerance: " << scientific << tol << "\n";

        cout << "Max iterations: " << maxIter << "\n";
        fout << "Max iterations: " << maxIter << "\n";

        set<double> roots;
        int count = 1;

        for(double x = start; x <= end; x += 0.5)
        {
            double curr = x, next;
            int iter;

            for(iter = 1; iter <= maxIter; iter++)
            {
                double fx = polyVal(coef, curr);
                double fdx = polyDeriv(coef, curr);

                if(fabs(fdx) < 1e-12) break;

                next = curr - fx / fdx;

                if(fabs(next - curr) < tol || fabs(fx) < tol) break;

                curr = next;
            }

            double r = round(next * 1000.0) / 1000.0;
            if(!isfinite(r)) continue;

            if(roots.find(r) == roots.end())
            {
                roots.insert(r);
                cout << "Root " << count << " = " << fixed << setprecision(6) << r
                     << " (Iterations: " << iter << ")\n";
                fout << "Root " << count << " = " << fixed << setprecision(6) << r
                     << " (Iterations: " << iter << ")\n";
                count++;
            }
        }

        if(roots.empty())
        {
            cout << "No real roots found.\n";
            fout << "No real roots found.\n";
        }

        cout << "\n";
        fout << "\n";
    }

    fin.close();
    fout.close();
    return 0;
}

```

#### Newton Raphson Input
```
2
3
1 -6 11 -6
0 3 0.001
50
2
1 -4 4
0 3 0.001
50

```

#### Newton Raphson Output
```
TestCase 1
Polynomial: f(x)=1x^3 - 6x^2 + 11x - 6
Search interval: [0, 3]
Tolerance: 1.000000e-03
Max iterations: 50
Root 1 = 1.000000 (Iterations: 5)
Root 2 = 2.000000 (Iterations: 5)
Root 3 = 3.000000 (Iterations: 5)

TestCase 2
Polynomial: f(x)=1x^2 - 4x + 4
Search interval: [0, 3]
Tolerance: 1.000000e-03
Max iterations: 50
Root 1 = 2.000000 (Iterations: 6)

```
[Back to Table of Contents](#table-of-contents)

---

## Solution of Linear Equations

### Gauss Elimination Method

#### Gauss Elimination Theory
The Gauss Elimination Method is used to solve a system of linear equations by transforming the system into an equivalent upper triangular form using matrix operations. 

A system of linear equations

a11 x1 + a12 x2 + … + a1n xn = b1

a21 x1 + a22 x2 + … + a2n xn = b2

.
.

an1 x1 + an2 x2 + … + ann xn = bn

is written in augmented matrix form as

[ a11 a12 … a1n | b1 ]

[ a21 a22 … a2n | b2 ]

[ .   .   .   . | .  ]

[ an1 an2 … ann | bn ] 

By applying elementary row operations, the augmented matrix is reduced to an upper triangular (row-echelon) form

[ a11  a12  …  a1n  | b1  ]

[ 0    a22' …  a2n' | b2' ]

[ 0    0    …  a3n' | b3' ]

[ 0    0    …  ann' | bn' ] 

Once the matrix is in upper triangular form, the solution of the system is obtained by back substitution, starting from the last equation and proceeding upward. 


#### Gauss Elimination Code
```cpp
#include <bits/stdc++.h>
using namespace std;

const double EPS = 1e-9;

int main()
{
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    int T;
    fin>>T;

    for(int i = 1; i <= T; i++)
    {
    int n;
    fin >> n;

    vector<vector<double>> a(n, vector<double>(n + 1));
    for (int i = 0; i < n; i++)
        for (int j = 0; j <= n; j++) fin >> a[i][j];
    
    fout << "Test Case " << i << "(" << n << "X" << n << ")" << endl;
    fout << fixed << setprecision(3);

    fout<<"\nMatrix:\n";
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j <= n; j++) fout << a[i][j] << "\t";
        fout<<endl;
    }
      
    for (int i = 0; i < n; i++)
    {
        int pivot = i;
        for (int k = i + 1; k < n; k++)
            if (fabs(a[k][i]) > fabs(a[pivot][i])) pivot = k;

        swap(a[i], a[pivot]);

        if (fabs(a[i][i]) < EPS) continue;

        for (int j = i + 1; j < n; j++)
        {
            double factor = a[j][i] / a[i][i];
            for (int k = i; k <= n; k++) a[j][k] -= factor * a[i][k];
        }
    }

    int rank = 0;
    bool noSolution = false;
    for (int i = 0; i < n; i++)
    {
        bool allZero = true;
        for (int j = 0; j < n; j++)
            if (fabs(a[i][j]) > EPS) allZero = false;

        if (allZero && fabs(a[i][n]) > EPS) noSolution = true;
        if (!allZero) rank++;
    }

    if (noSolution) fout << "\nThe system has NO solution.\n";
    else if (rank < n) fout << "\nThe system has INFINITE solutions.\n";
    else
    {
        vector<double> x(n);
        for (int i = n - 1; i >= 0; i--)
        {
            x[i] = a[i][n];
            for (int j = i + 1; j < n; j++) x[i] -= a[i][j] * x[j];
            x[i] /= a[i][i];
        }

        fout << "\nThe system has a UNIQUE solution:\n";
        for (int i = 0; i < n; i++) fout << "x" << i + 1 << " = " << x[i] << "\n";
    }

    if(i < T) fout << endl;
    }

    fin.close();
    fout.close();

    return 0;
}

```

#### Gauss Elimination Input
```
4
3
2 6 0 -11
6 20 -6 -3
0 6 -18 -1

3
5 3 7 4
3 26 2 9
7 2 10 5

3
2 3 4 11
1 5 7 15
3 11 13 25

2
1 1 2
1 1 3
```

#### Gauss Elimination Output
```
Test Case 1(3X3)

Matrix:
2.000	6.000	0.000	-11.000	
6.000	20.000	-6.000	-3.000	
0.000	6.000	-18.000	-1.000	

The system has NO solution.

Test Case 2(3X3)

Matrix:
5.000	3.000	7.000	4.000	
3.000	26.000	2.000	9.000	
7.000	2.000	10.000	5.000	

The system has INFINITE solutions.

Test Case 3(3X3)

Matrix:
2.000	3.000	4.000	11.000	
1.000	5.000	7.000	15.000	
3.000	11.000	13.000	25.000	

The system has a UNIQUE solution:
x1 = 2.000
x2 = -3.000
x3 = 4.000

Test Case 4(2X2)

Matrix:
1.000	1.000	2.000	
1.000	1.000	3.000	

The system has NO solution.

```
[Back to Table of Contents](#table-of-contents)

---

### Gauss Jordan Elimination Method

#### Gauss Jordan Theory
The Gauss–Jordan Elimination method is a direct numerical technique used to solve a system of linear equations. It is an extension of the Gauss Elimination method and transforms the given system into a reduced row-echelon form (RREF), where the solution can be obtained directly.

A system of linear equations can be represented in matrix form as:

A X = B

where
A is the coefficient matrix,
X is the column matrix of unknowns,
B is the constant column matrix.

This system is written as an augmented matrix:

[A | B]

The objective of the Gauss–Jordan method is to convert the augmented matrix into the form:

[I | X]

where I is the identity matrix and X contains the solutions of the variables.

This transformation is achieved using elementary row operations:

Interchanging two rows

Multiplying a row by a non-zero constant

Adding or subtracting a multiple of one row to another row

In Gauss–Jordan elimination, each pivot element (leading diagonal element) is first converted to 1, and then all other elements in that column are made zero. This process is repeated for every variable until the left side of the augmented matrix becomes an identity matrix.

Unlike Gauss Elimination, which produces an upper triangular matrix and requires back substitution, Gauss–Jordan elimination directly yields the solution without back substitution.

The method is widely used for:

Solving systems of linear equations

Finding the inverse of a matrix

Checking consistency of linear systems

Although Gauss–Jordan elimination requires more computations than Gauss elimination, it provides a straightforward and systematic approach to obtaining exact solutions.

#### Gauss Jordan Code
```cpp
#include <bits/stdc++.h>
using namespace std;

const double EPS = 1e-9;

int main()
{  
    ifstream fin("input.txt");
    ofstream fout("output.txt");
    
    int T;
    fin>>T;

    for(int i = 1; i <= T; i++)
    {
    int n;
    fin >> n;

    vector<vector<double>> a(n, vector<double>(n + 1));

    for (int i = 0; i < n; i++)
        for (int j = 0; j <= n; j++) fin >> a[i][j];

    fout<<"Test Case "<<i<<"("<<n<<"X"<<n<<")"<<"\n\n";
    fout<<fixed<<setprecision(3);

    fout<<"Matrix:\n";
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j <= n; j++) fout << a[i][j] << "\t";
        fout << endl;
    }

    int rank = 0;

    for (int col = 0; col < n; col++)
    {
        int pivot = rank;
        for (int i = rank; i < n; i++)
        {
            if (fabs(a[i][col]) > fabs(a[pivot][col])) pivot = i;
        }

        if (fabs(a[pivot][col]) < EPS) continue;

        swap(a[pivot], a[rank]);

        double div = a[rank][col];
        for (int j = 0; j <= n; j++) a[rank][j] /= div;

        for (int i = 0; i < n; i++)
        {
            if (i != rank)
            {
                double factor = a[i][col];
                for (int j = 0; j <= n; j++) a[i][j] -= factor * a[rank][j];
            }
        }
        rank++;
    }

    bool noSolution = false;
    for (int i = 0; i < n; i++)
    {
        bool allZero = true;
        for (int j = 0; j < n; j++)
            if (fabs(a[i][j]) > EPS) allZero = false;
            if (allZero && fabs(a[i][n]) > EPS) noSolution = true;
    }

    fout<<endl;

    fout<<"Reduced Row Echelon Form:\n";
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j <= n; j++) fout << a[i][j]<<"\t";
        fout<<endl;
    }

    if (noSolution) fout <<"\nThe system has NO solution(Inconsistent).\n";
    else if (rank < n) fout <<"\nThe system has INFINITE solutions.\n";
    else
    {
        fout << "\nSolution:\n";
        for (int i = 0; i < n; i++) fout << "x" << i + 1 << " = "<< a[i][n] << "\n";
        fout << "\nThe system has a UNIQUE solution.\n";
    }

    if(i < T) fout << endl;
    }

    fin.close();
    fout.close();

    return 0;
}

```

#### Gauss Jordan Input
```
4
3
2 6 0 -11
6 20 -6 -3
0 6 -18 -1

3
5 3 7 4
3 26 2 9
7 2 10 5

3
2 3 4 11
1 5 7 15
3 11 13 25

2
1 1 2
1 1 3
```

#### Gauss Jordan Output
```
Test Case 1(3X3)

Matrix:
2.000	6.000	0.000	-11.000	
6.000	20.000	-6.000	-3.000	
0.000	6.000	-18.000	-1.000	

Reduced Row Echelon Form:
1.000	0.000	9.000	0.056	
0.000	1.000	-3.000	-0.167	
0.000	0.000	-0.000	-10.111	

The system has NO solution(Inconsistent).

Test Case 2(3X3)

Matrix:
5.000	3.000	7.000	4.000	
3.000	26.000	2.000	9.000	
7.000	2.000	10.000	5.000	

Reduced Row Echelon Form:
1.000	0.000	1.455	0.636	
0.000	1.000	-0.091	0.273	
0.000	0.000	-0.000	-0.000	

The system has INFINITE solutions.

Test Case 3(3X3)

Matrix:
2.000	3.000	4.000	11.000	
1.000	5.000	7.000	15.000	
3.000	11.000	13.000	25.000	

Reduced Row Echelon Form:
1.000	0.000	0.000	2.000	
-0.000	1.000	0.000	-3.000	
0.000	0.000	1.000	4.000	

Solution:
x1 = 2.000
x2 = -3.000
x3 = 4.000

The system has a UNIQUE solution.

Test Case 4(2X2)

Matrix:
1.000	1.000	2.000	
1.000	1.000	3.000	

Reduced Row Echelon Form:
1.000	1.000	2.000	
0.000	0.000	1.000	

The system has NO solution(Inconsistent).

```
[Back to Table of Contents](#table-of-contents)

---

### LU Decomposition Method

#### LU Decomposition Theory
LU Decomposition is a numerical method used to solve a system of linear equations by factorizing the coefficient matrix into a lower triangular matrix and an upper triangular matrix. It simplifies the solution process by converting the original system into two triangular systems.

Matrix Representation

Given the system:

A X = B

The coefficient matrix A is decomposed as:

A = L U

where
L = lower triangular matrix (with diagonal elements equal to 1)
U = upper triangular matrix

Thus,

L U X = B

Solve in two steps:

L Y = B
U X = Y

Algorithm (Short Steps)

Read the coefficient matrix A and constant vector B.

Decompose matrix A into L and U such that A = L U.

Use forward substitution to solve L Y = B.

Use backward substitution to solve U X = Y.

Obtain the solution vector X.

If you want the same format for Gauss Elimination / Gauss–Jordan / Cholesky, just tell me 👍

#### LU Decomposition Code
```cpp
#include <bits/stdc++.h>
using namespace std;

const double EPS = 1e-9;

void LUdecomposition(vector<vector<double>> A, vector<vector<double>> &L, vector<vector<double>> &U, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int k = i; k < n; k++)
        {
            double sum = 0;
            for (int j = 0; j < i; j++)
                sum += (L[i][j] * U[j][k]);
            U[i][k] = A[i][k] - sum;
        }

        for (int k = i; k < n; k++)
        {
            if (i == k)
                L[i][i] = 1;
            else {
                double sum = 0;
                for (int j = 0; j < i; j++)
                    sum += (L[k][j] * U[j][i]);
                if (fabs(U[i][i]) < 1e-9) cout << "Division by zero detected during LU decomposition!" << endl;
                L[k][i] = (A[k][i] - sum) / U[i][i];
            }
        }
    }
}

vector<double> forwardSub(const vector<vector<double>>& L, const vector<double>& b)
{
    int n = b.size();
    vector<double> y(n);
    for (int i = 0; i < n; i++)
    {
        double sum = 0;
        for (int j = 0; j < i; j++) sum += L[i][j] * y[j];
        y[i] = b[i] - sum;
    }
    return y;
}

vector<double> backwardSub(const vector<vector<double>>& U, const vector<double>& y)
{
    int n = y.size();
    vector<double> x(n);
    for (int i = n - 1; i >= 0; i--)
    {
        double sum = 0;
        for (int j = i + 1; j < n; j++) sum += U[i][j] * x[j];
        x[i] = (y[i] - sum) / U[i][i];
    }
    return x;
}

int matrixRank(vector<vector<double>> M)
{
    if (M.empty()) return 0;
    int rows = M.size();
    int cols = M[0].size();
    int rank = 0;

    for (int col = 0; col < cols; col++)
    {
        int pivot = -1;
        for (int i = rank; i < rows; i++)
        {
            if (col < M[i].size() && fabs(M[i][col]) > EPS)
            {
                pivot = i;
                break;
            }
        }
        if (pivot == -1) continue;

        swap(M[pivot], M[rank]);
        double div = M[rank][col];
        for (int j = col; j < cols; j++)
            M[rank][j] /= div;

        for (int i = 0; i < rows; i++)
        {
            if (i != rank)
            {
                double factor = (col < M[i].size()) ? M[i][col] : 0;
                for (int j = col; j < cols; j++)
                    if (j < M[i].size()) M[i][j] -= factor * M[rank][j];
            }
        }
        rank++;
    }
    return rank;
}

int main()
{
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    int T;
    fin>>T;

    for(int i = 1; i <= T; i++)
    {
    int n;
    fin>>n;

    vector<vector<double>> A(n, vector<double>(n));
    vector<double> B(n);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++) fin >> A[i][j];
        fin >> B[i];
    }

    vector<vector<double>> Aug(n, vector<double>(n+1));
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++) Aug[i][j] = A[i][j];
        Aug[i][n] = B[i];
    }

    vector<vector<double>> L(n, vector<double>(n, 0));
    vector<vector<double>> U(n, vector<double>(n, 0));

    int rankA = matrixRank(A);
    int rankAug = matrixRank(Aug);

    fout << "Test Case " << i << endl;
    fout << fixed << setprecision(3); 

    fout<<"Augmented Matrix:\n";
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++) fout << A[i][j] << "\t";
        fout << B[i] <<endl;
    }

    LUdecomposition(A, L, U, n);

    fout << "\nLower Triangular Matrix (L):\n";
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++) fout << L[i][j] << "\t";
        fout << endl;
    }

    fout << "\nUpper Triangular Matrix (U):\n";
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++) fout << U[i][j] << "\t";
        fout << endl;
    }


    if (rankA != rankAug) fout << "\nThe system is INCONSISTENT (No solution).\n";
    else if (rankA < n) fout << "\nThe system has INFINITE solutions.\n";
    else
    {
        vector<double> y = forwardSub(L, B);
        vector<double> x = backwardSub(U, y);  

        fout<< "\nValues of Y:\n";
        for(int i = 0; i < n; i++) fout << "y" << i+1 << " = "<< y[i] <<endl;

        fout << "\nUnique solution:\n";
        for (int i = 0; i < n; i++) fout << "x" << i+1 << " = " << x[i] << endl;
    }

    if(i < T) fout<<endl;
    }

    fin.close();
    fout.close();

    return 0;
}
```

#### LU Decomposition Input
```
5
3
2 6 0 -11
6 20 -6 -3
0 6 -18 -1

3
5 3 7 4
3 26 2 9
7 2 10 5

3
2 3 4 11
1 5 7 15
3 11 13 25

2
1 1 2
1 1 3

3
2 1 -1 8
-3 -1 2 -11
-2 1 2 -3
```

#### LU Decomposition Output
```
Test Case 1
Augmented Matrix:
2.000	6.000	0.000	-11.000
6.000	20.000	-6.000	-3.000
0.000	6.000	-18.000	-1.000

Lower Triangular Matrix (L):
1.000	0.000	0.000	
3.000	1.000	0.000	
0.000	3.000	1.000	

Upper Triangular Matrix (U):
2.000	6.000	0.000	
0.000	2.000	-6.000	
0.000	0.000	0.000	

The system is INCONSISTENT (No solution).

Test Case 2
Augmented Matrix:
5.000	3.000	7.000	4.000
3.000	26.000	2.000	9.000
7.000	2.000	10.000	5.000

Lower Triangular Matrix (L):
1.000	0.000	0.000	
0.600	1.000	0.000	
1.400	-0.091	1.000	

Upper Triangular Matrix (U):
5.000	3.000	7.000	
0.000	24.200	-2.200	
0.000	0.000	0.000	

The system has INFINITE solutions.

Test Case 3
Augmented Matrix:
2.000	3.000	4.000	11.000
1.000	5.000	7.000	15.000
3.000	11.000	13.000	25.000

Lower Triangular Matrix (L):
1.000	0.000	0.000	
0.500	1.000	0.000	
1.500	1.857	1.000	

Upper Triangular Matrix (U):
2.000	3.000	4.000	
0.000	3.500	5.000	
0.000	0.000	-2.286	

Values of Y:
y1 = 11.000
y2 = 9.500
y3 = -9.143

Unique solution:
x1 = 2.000
x2 = -3.000
x3 = 4.000

Test Case 4
Augmented Matrix:
1.000	1.000	2.000
1.000	1.000	3.000

Lower Triangular Matrix (L):
1.000	0.000	
1.000	1.000	

Upper Triangular Matrix (U):
1.000	1.000	
0.000	0.000	

The system is INCONSISTENT (No solution).

Test Case 5
Augmented Matrix:
2.000	1.000	-1.000	8.000
-3.000	-1.000	2.000	-11.000
-2.000	1.000	2.000	-3.000

Lower Triangular Matrix (L):
1.000	0.000	0.000	
-1.500	1.000	0.000	
-1.000	4.000	1.000	

Upper Triangular Matrix (U):
2.000	1.000	-1.000	
0.000	0.500	0.500	
0.000	0.000	-1.000	

Values of Y:
y1 = 8.000
y2 = 1.000
y3 = 1.000

Unique solution:
x1 = 2.000
x2 = 3.000
x3 = -1.000
```
[Back to Table of Contents](#table-of-contents)

---

### Matrix Inversion Method

#### Matrix Inversion Theory
The Matrix Inversion Method is a numerical technique used to solve a system of linear equations by finding the inverse of the coefficient matrix. Once the inverse matrix is obtained, the solution can be calculated directly using matrix multiplication.

Matrix Representation

Given the system:

A X = B

If the inverse of A exists, then:

A⁻¹ A X = A⁻¹ B

⇒ X = A⁻¹ B

where
A⁻¹ is the inverse of the coefficient matrix A.

Algorithm (Short Steps):

1.Read the coefficient matrix A and constant vector B.

2.Check whether det(A) ≠ 0 (inverse exists).

3.Find the inverse of matrix A using Gauss–Jordan elimination.

4.Multiply A⁻¹ with B to obtain X.

5.The resulting vector X is the solution of the system.

#### Matrix Inversion Code
```cpp
#include <bits/stdc++.h>
using namespace std;

const double EPS = 1e-9;

void printMatrix(ofstream& fout,vector<vector<double>> &A)
{
    for (auto &row : A)
    {
        for (double v : row)
        fout << setw(10) << fixed << setprecision(3) << v << " ";
        fout <<endl;
    }
    fout <<endl;
}

double Determinant(vector<vector<double>> A, int n)
{
    double det = 1;

    for (int i = 0; i < n; i++)
    {
        int pivot = i;
        for (int j = i + 1; j < n; j++)
        if (fabs(A[j][i]) > fabs(A[pivot][i])) pivot = j;

        if (fabs(A[pivot][i]) < EPS) return 0;

        if (pivot != i)
        {
            swap(A[pivot], A[i]);
            det *= -1;
        }

        det *= A[i][i];

        for (int j = i + 1; j < n; j++)
        {
            double factor = A[j][i] / A[i][i];
            for (int k = i; k < n; k++) A[j][k] -= factor * A[i][k];
        }
    }
    return det;
}

void CheckConsistency(ofstream& fout,vector<vector<double>> A, int n)
{
    int m = n + 1;

    for (int i = 0; i < n; i++)
    {
        int pivot = i;
        for (int j = i; j < n; j++)
        if (fabs(A[j][i]) > fabs(A[pivot][i])) pivot = j;

        if (fabs(A[pivot][i]) < EPS) continue;

        swap(A[pivot], A[i]);

        for (int j = i + 1; j < n; j++)
        {
            double factor = A[j][i] / A[i][i];
            for (int k = i; k < m; k++) A[j][k] -= factor * A[i][k];
        }
    }

    for (int i = 0; i < n; i++)
    {
        bool allZero = true;
        for (int j = 0; j < n; j++)
            if (fabs(A[i][j]) > EPS) allZero = false;

        if (allZero && fabs(A[i][n]) > EPS)
        {
            fout <<"No solution (Inconsistent system)"<<endl;
            return;
        }
    }
    fout <<"Infinitely many solutions (Consistent system)"<<endl;
}

vector<vector<double>> Transpose(vector<vector<double>> &A, int n) {
    vector<vector<double>> T(n, vector<double>(n));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) T[j][i] = A[i][j];
    return T;
}

vector<vector<double>> Inverse(vector<vector<double>> &A, int n) {
    vector<vector<double>> I(n, vector<double>(2 * n));

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++) I[i][j] = A[i][j];
        I[i][i + n] = 1;
    }

    for (int i = 0; i < n; i++)
    {
        int pivot = i;
        for (int j = i; j < n; j++)
            if (fabs(I[j][i]) > fabs(I[pivot][i])) pivot = j;

        swap(I[pivot], I[i]);

        double div = I[i][i];
        for (int j = 0; j < 2 * n; j++) I[i][j] /= div;

        for (int j = 0; j < n; j++)
        {
            if (j == i) continue;
            double factor = I[j][i];
            for (int k = 0; k < 2 * n; k++) I[j][k] -= factor * I[i][k];
        }
    }

    vector<vector<double>> inv(n, vector<double>(n));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) inv[i][j] = I[i][j + n];

    return inv;
}

vector<double> SolveX(vector<vector<double>> &inv, vector<double> &B, int n)
{
    vector<double> X(n, 0);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) X[i] += inv[i][j] * B[j];
    return X;
}

int main()
{
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    int T;
    fin>>T;

    for(int i = 1; i <= T; i++)
    {
    int n;
    fin >> n;

    vector<vector<double>> A(n, vector<double>(n));
    vector<vector<double>> Aug(n, vector<double>(n + 1));
    vector<double> B(n);

    for (int i = 0; i < n; i++)
        for (int j = 0; j <= n; j++) fin >> Aug[i][j];

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++) A[i][j] = Aug[i][j];
        B[i] = Aug[i][n];
    }

    fout << "Test Case " << i <<endl;

    double det = Determinant(A, n);
    fout << "\nDeterminant = " << det << "\n\n";

    if (fabs(det) < EPS)
    {
        fout << "Determinant is zero\n";
        CheckConsistency(fout, Aug, n);
    }
    else
    {
        fout << "Transpose of A:\n";
        auto T = Transpose(A, n);
        printMatrix(fout, T);

        fout << "Inverse of A:\n";
        auto I = Inverse(A, n);
        printMatrix(fout, I);

        fout << "Solution Vector X:\n";
        auto X = SolveX(I, B, n);
        for (int i = 0; i < n; i++)
            fout << "x" << i + 1 << " = " << X[i] <<endl;
    }

    if(i < T) fout<<endl;
    }

    fin.close();
    fout.close();

    return 0;
}
```

#### Matrix Inversion Input
```
3
3
2 6 0 -11
6 20 -6 -3
0 6 -18 -1

3
5 3 7 4
3 26 2 9
7 2 10 5

3
2 3 4 11
1 5 7 15
3 11 13 25
```

#### Matrix Inversion Output
```
Test Case 1

Determinant = 0

Determinant is zero
No solution (Inconsistent system)

Test Case 2

Determinant = 0

Determinant is zero
Infinitely many solutions (Consistent system)

Test Case 3

Determinant = -16

Transpose of A:
     2.000      1.000      3.000 
     3.000      5.000     11.000 
     4.000      7.000     13.000 

Inverse of A:
     0.750     -0.312     -0.063 
    -0.500     -0.875      0.625 
     0.250      0.812     -0.437 

Solution Vector X:
x1 = 2.000
x2 = -3.000
x3 = 4.000

```
[Back to Table of Contents](#table-of-contents)

---

## Solution of Differential Equations

### Runge-Kutta 4th Order Method

#### Runge-Kutta 4th Order Theory


The Runge-Kutta 4th Order (RK4) method is a numerical technique used to
solve first-order ordinary differential equations of the form:

dy/dx = f(x, y)

Given an initial value y(x₀) = y₀, RK4 estimates the value of y at the
next point xₙ₊₁ = xₙ + h using four intermediate slopes:

k₁ = h · f(xₙ, yₙ)  
k₂ = h · f(xₙ + h/2, yₙ + k₁/2)  
k₃ = h · f(xₙ + h/2, yₙ + k₂/2)  
k₄ = h · f(xₙ + h, yₙ + k₃)

The next value of y is calculated as:

yₙ₊₁ = yₙ + (k₁ + 2k₂ + 2k₃ + k₄) / 6

RK4 provides high accuracy and is widely used in engineering and
scientific computations.


#### Runge-Kutta 4th Order Code
```cpp
#include <bits/stdc++.h>
using namespace std;

float dydx(float x, float y)
{
    return (x - y) / 2;
}


float rungeKutta(float x0, float y0, float x, float h)
{
    int n = (int)((x - x0) / h);

    float k1, k2, k3, k4;
    float y = y0;

    for (int i = 1; i <= n; i++)
    {
        k1 = h * dydx(x0, y);
        k2 = h * dydx(x0 + h / 2, y + k1 / 2);
        k3 = h * dydx(x0 + h / 2, y + k2 / 2);
        k4 = h * dydx(x0 + h, y + k3);

        y = y + (k1 + 2*k2 + 2*k3 + k4) / 6;
        x0 = x0 + h;
    }

    return y;
}

int main()
{
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    if (!fin)
    {
        cout << "Error: input.txt file not found!" << endl;
        return 1;
    }

    float x0, y0, x, h;
    fin >> x0 >> y0 >> x >> h;

    float result = rungeKutta(x0, y0, x, h);
    fout << fixed << setprecision(4);
    fout << "Initial value of x (x0)   : " << x0 << endl;
    fout << "Initial value of y (y0)   : " << y0 << endl;
    fout << "Value of x required       : " << x  << endl;
    fout << "Step size (h)             : " << h  << endl;

    
    fout << " Result\n";
    fout << "Approximate value of y at x = " << x
         << " is: " << result << endl;
    

    fin.close();
    fout.close();

    return 0;
}


```

#### Runge-Kutta 4th Order Input
```
0 1 2 0.2
```


#### Runge-Kutta 4th Order Output
```
Initial value of x (x0)   : 0.0000
Initial value of y (y0)   : 1.0000
Value of x required       : 2.0000
Step size (h)             : 0.2000
 Result
Approximate value of y at x = 2.0000 is: 1.1036

```
[Back to Table of Contents](#table-of-contents)

---

## Numerical Integration

### Simpson's 1/3 Rule

#### Simpson's 1/3 Theory

Simpson’s 1/3 Rule is a numerical integration technique used to approximate the definite integral of a function. It provides a more accurate result than the Trapezoidal Rule by approximating the integrand with a second-degree polynomial (parabola) over each pair of subintervals. This method is particularly effective for smooth functions and is exact for polynomials of degree three or less.

Mathematical Principle
Let f(x) be a continuous function on the interval [a, b]. The interval is divided into an even number (n) of equal subintervals, each of width:
h = (b - a) / n

The approximate value of the integral is given by:
∫[a to b] f(x) dx ≈ (h / 3) [ f(x0) + f(xn)
+ 4 (f(x1) + f(x3) + ... + f(xn-1))
+ 2 (f(x2) + f(x4) + ... + f(xn-2)) ]

where n must be even.

Algorithm (Simpson’s 1/3 Rule)
1. Choose the limits of integration a and b.
2. Select an even number of subintervals n.
3. Compute the step size h = (b - a) / n.
4. Evaluate the function at each point x_i = a + i*h.
5. Apply Simpson’s 1/3 formula to calculate the approximate integral.
6. The result obtained is the numerical approximation of the definite integral.

Features
- Higher accuracy than Trapezoidal Rule
- Exact for cubic polynomials
- Requires even number of subintervals
- Suitable for smooth functions
- Simple and efficient implementation
- Default error tolerance depends on step size


#### Simpson's 1/3 Code
```cpp
#include <bits/stdc++.h>
using namespace std;

using db = double;

db evalPoly(db x, const vector<db>& coeff)
{
    db sum = 0;
    db powx = 1;
    for(db c : coeff)
    {
        sum += c * powx;
        powx *= x;
    }
    return sum;
}

void simpsonRule(int tc, ifstream &fin, ofstream &fout)
{
    int deg;
    fin >> deg;
    vector<db> coef(deg+1);
    for(int i=0; i<=deg; i++) fin >> coef[i];

    db a, b;
    int n;
    fin >> a >> b >> n;

    if(n % 2 != 0)
    {
        cout << "Test Case #" << tc << ": Error - n must be even!\n\n";
        fout << "Test Case #" << tc << ": Error - n must be even!\n\n";
        return;
    }

    db h = (b - a) / n;
    db sum = evalPoly(a, coef) + evalPoly(b, coef);

    for(int i=1; i<n; i++)
    {
        db x = a + i*h;
        sum += (i % 2 == 1) ? 4*evalPoly(x, coef) : 2*evalPoly(x, coef);
    }

    db result = (h/3) * sum;

    cout << "Test Case #" << tc << "\n";
    cout << "Degree: " << deg << "\n";
    cout << "Coefficients: ";
    for(db c : coef) cout << c << " ";
    cout << "\nInterval: [" << a << ", " << b << "]\n";
    cout << "Intervals (n): " << n << "\n";
    cout << "Step size (h): " << h << "\n";
    cout << "Integral result: " << result << "\n\n";

    fout << "Test Case #" << tc << "\n";
    fout << "Degree: " << deg << "\n";
    fout << "Coefficients: ";
    for(db c : coef) fout << c << " ";
    fout << "\nInterval: [" << a << ", " << b << "]\n";
    fout << "Intervals (n): " << n << "\n";
    fout << "Step size (h): " << h << "\n";
    fout << "Integral result: " << result << "\n\n";
}

int main()
{
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    if(!fin)
    {
        cout << "Error: input.txt not found!\n";
        return 0;
    }

    int t;
    fin >> t;
    cout << "Total Test Cases: " << t << "\n\n";
    fout << "Total Test Cases: " << t << "\n\n";

    for(int tc = 1; tc <= t; tc++)
        simpsonRule(tc, fin, fout);

    cout << "All results written to output.txt!\n";
    return 0;
}

```

#### Simpson's 1/3 Input
```
2
3
1 0 -2 1
0 2 4
2
2
2 -3 1
-1 1 6

```

#### Simpson's 1/3 Output
```
Test Case #1
Degree: 3
Coefficients: 1 0 -2 1 
Interval: [0, 2]
Intervals (n): 4
Step size (h): 0.5
Integral result: 0.666667

Test Case #2
Degree: 2
Coefficients: 2 -3 1 
Interval: [-1, 1]
Intervals (n): 6
Step size (h): 0.333333
Integral result: 0.000000

```
[Back to Table of Contents](#table-of-contents)

---

### Simpson's 3/8 Rule

#### Simpson's 3/8 Theory
Simpson’s 3/8 Rule
Theory

Simpson’s 3/8 Rule is a numerical integration technique used to approximate the definite integral of a function. It is an extension of Simpson’s 1/3 Rule and approximates the integrand using a third-degree polynomial (cubic) over three consecutive subintervals. This method provides good accuracy for smooth functions and is exact for polynomials of degree three or less.

Algorithm (Simpson’s 3/8 Rule)
1. Choose the limits of integration a and b.
2. Select the number of subintervals n such that n is a multiple of 3.
3. Compute the step size h = (b-a)/h.
4. Evaluate the function at each point xi = a + i*h.
5. Apply Simpson’s 3/8 formula to calculate the approximate integral.
6. The result obtained is the numerical approximation of the definite integral.

#### Simpson's 3/8 Code
```python
#include <bits/stdc++.h>
using namespace std;

using db = double;

db evalPoly(db x, const vector<db>& coeff)
{
    db sum = 0, powx = 1;
    for(db c : coeff)
    {
        sum += c * powx;
        powx *= x;
    }
    return sum;
}

void simpson38(int tc, ifstream &fin, ofstream &fout)
{
    int deg;
    fin >> deg;
    vector<db> coef(deg + 1);
    for(int i=0; i<=deg; i++) fin >> coef[i];

    db a, b;
    int n;
    fin >> a >> b >> n;

    if(n % 3 != 0)
    {
        cout << "Test Case #" << tc << ": Error - n must be a multiple of 3!\n\n";
        fout << "Test Case #" << tc << ": Error - n must be a multiple of 3!\n\n";
        return;
    }

    db h = (b - a) / n;
    db sum = evalPoly(a, coef) + evalPoly(b, coef);

    for(int i=1; i<n; i++)
    {
        db x = a + i*h;
        sum += (i % 3 == 0) ? 2*evalPoly(x, coef) : 3*evalPoly(x, coef);
    }

    db result = (3*h/8) * sum;

    cout << "Test Case #" << tc << "\n";
    cout << "Degree: " << deg << "\n";
    cout << "Coefficients: ";
    for(db c : coef) cout << c << " ";
    cout << "\nInterval: [" << a << ", " << b << "]\n";
    cout << "Intervals (n): " << n << "\n";
    cout << "Step size (h): " << h << "\n";
    cout << "Integral result: " << result << "\n\n";

    fout << "Test Case #" << tc << "\n";
    fout << "Degree: " << deg << "\n";
    fout << "Coefficients: ";
    for(db c : coef) fout << c << " ";
    fout << "\nInterval: [" << a << ", " << b << "]\n";
    fout << "Intervals (n): " << n << "\n";
    fout << "Step size (h): " << h << "\n";
    fout << "Integral result: " << result << "\n\n";
}

int main()
{
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    if(!fin)
    {
        cout << "Error: input.txt not found!\n";
        return 0;
    }

    int t;
    fin >> t;
    cout << "Total Test Cases: " << t << "\n\n";
    fout << "Total Test Cases: " << t << "\n\n";

    for(int tc=1; tc<=t; tc++)
        simpson38(tc, fin, fout);

    cout << "All results written to output.txt!\n";
    return 0;
}

```

#### Simpson's 3/8 Input
```
2
3
1 0 -2 1
0 2 6
2
2
2 -3 1
-1 1 3

```

#### Simpson's 3/8 Output
```
Test Case #1
Degree: 3
Coefficients: 1 0 -2 1 
Interval: [0, 2]
Intervals (n): 6
Step size (h): 0.333333
Integral result: 0.666667

Test Case #2
Degree: 2
Coefficients: 2 -3 1 
Interval: [-1, 1]
Intervals (n): 3
Step size (h): 0.666667
Integral result: 0.000000

```
[Back to Table of Contents](#table-of-contents)

---

## Interpolation Methods

### Newton Forward Interpolation Method
#### Newton Forward Theory
Newton's forward interpolation is a numerical method used to approximate a function value near the beginning of a dataset. It strictly requires the data points to be spaced at equal intervals.​

The Formula
The method uses "forward differences" (
Δ
Δ) derived from the dataset to build a polynomial. The value 
y
y at a specific 
x
x is:

y(x)=y0+pΔy0+((p(p−1))/2!) *Δ^ 2y0 + ((p(p−1)(p−2))/3!)* Δ^3y0+…

Short Algorithm: Newton’s Forward Interpolation

1. Initialize: Read arrays X and Y, and the target value x.
2. Calculate p: Compute step size h = X[1] - X[0] and p = (x - X[0]) / h.
3. Build Table: Construct a difference table where every new column is the difference of the previous one: diff[i][j] = diff[i+1][j-1] - diff[i][j-1].
4. Compute Sum: Initialize sum = Y[0] and u = 1. Loop from j = 1 to n-1:
   - Update u = u * (p - (j - 1)).
   - Add term: sum = sum + (u * diff[0][j]) / factorial(j).
5. Output: Return sum.
``

#### Newton Forward Code
```cpp
#include <bits/stdc++.h>
using namespace std;


float u_cal(float u, int n)
{
    float res = u;
    for (int i = 1; i < n; i++)
        res *= (u - i);
    return res;
}


int fact(int n)
{
    int f = 1;
    for (int i = 2; i <= n; i++)
        f *= i;
    return f;
}

int main()
{
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    if (!fin) {
        fout << "Error opening input file!" << endl;
        return 0;
    }

    int n;
    fin >> n;

    vector<float> x(n);
    for (int i = 0; i < n; i++)
        fin >> x[i];

    vector<vector<float>> y(n, vector<float>(n));

    for (int i = 0; i < n; i++)
        fin >> y[i][0];

    float value;
    fin >> value;

    fin.close();

    
    for (int i = 1; i < n; i++) {
        for (int j = 0; j < n - i; j++) {
            y[j][i] = y[j + 1][i - 1] - y[j][i - 1];
        }
    }

    
    fout << "Forward Difference Table\n\n";
    for (int i = 0; i < n; i++) {
        fout << setw(6) << x[i] << "\t";
        for (int j = 0; j < n - i; j++)
            fout << setw(8) << y[i][j] << "\t";
        fout << endl;
    }

    
    float u = (value - x[0]) / (x[1] - x[0]);
    float sum = y[0][0];

    for (int i = 1; i < n; i++) {
        sum += (u_cal(u, i) * y[0][i]) / fact(i);
    }

    fout << "\nInterpolated value at x = " << value << " is " << sum << endl;

    fout.close();

    return 0;
}

```
#### Newton Forward Input
```
4
45 50 55 60
0.7071 0.7660 0.8192 0.8660
52

```

#### Newton Forward Output
```
Forward Difference Table

    45	  0.7071	  0.0589	-0.00569999	-0.000699997	
    50	   0.766	  0.0532	-0.00639999	
    55	  0.8192	  0.0468	
    60	   0.866	

Interpolated value at x = 52 is 0.788003
```
[Back to Table of Contents](#table-of-contents)

---

### Newton Backward Interpolation Method
#### Newton Backward Theory
 Newton's Backward Interpolation Method

Newton's backward interpolation is a numerical method designed to approximate values located near the END of a dataset. It strictly requires the data points to be spaced at equal intervals.

Formula:
y(x) = y_n + p * (nabla y_n) + (p(p+1)/2!) * (nabla^2 y_n) + ...

ALGORITHM: Newton’s Forward Interpolation
1.Read the number of data points n.
2.Read the equally spaced values of x and corresponding values of y.
3.Construct the forward difference table using the given y values.
4.Read the value of x for which interpolation is required.
5.Compute the step size h = x1 − x0.
6.Compute u = (x − x0) / h.
7.Initialize result = y0.
8.For i = 1 to n − 1:
Compute term = (u(u−1)(u−2)…)/(i!)
9.Multiply term with the i-th forward difference.
10.Add the term to result.
11.Output the interpolated value.

#### Newton Backward Code
```cpp
#include <bits/stdc++.h>
using namespace std;

float u_cal(float u, int n)
{
    float res = u;
    for (int i = 1; i < n; i++)
        res *= (u + i);
    return res;
}


int fact(int n)
{
    int f = 1;
    for (int i = 2; i <= n; i++)
        f *= i;
    return f;
}

int main()
{
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    if (!fin) {
        fout << "Error opening input file!" << endl;
        return 0;
    }

    int n;
    fin >> n;

    vector<float> x(n);
    for (int i = 0; i < n; i++)
        fin >> x[i];

    vector<vector<float>> y(n, vector<float>(n));

    for (int i = 0; i < n; i++)
        fin >> y[i][0];

    float value;
    fin >> value;

    fin.close();

    
    for (int i = 1; i < n; i++) {
        for (int j = n - 1; j >= i; j--) {
            y[j][i] = y[j][i - 1] - y[j - 1][i - 1];
        }
    }


    fout << "Backward Difference Table\n\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++)
            fout << setw(8) << y[i][j] << "\t";
        fout << endl;
    }

    
    float u = (value - x[n - 1]) / (x[1] - x[0]);
    float sum = y[n - 1][0];

    for (int i = 1; i < n; i++) {
        sum += (u_cal(u, i) * y[n - 1][i]) / fact(i);
    }

    fout << "\nInterpolated value at x = " << value
         << " is " << sum << endl;

    fout.close();

    return 0;
}

````
#### Newton Backward Input
```
5
1891 1901 1911 1921 1931
46 66 81 93 101
1925
```

#### Newton Backward Output
```
Backward Difference Table

      46	
      66	      20	
      81	      15	      -5	
      93	      12	      -3	       2	
     101	       8	      -4	      -1	      -3	

Interpolated value at x = 1925 is 96.8368
````
[Back to Table of Contents](#table-of-contents)

---
### Newton Divided Difference Interpolation Method
#### Newton Divided Difference Theory
Newton’s Divided Difference Table

Newton’s divided difference method is an interpolation technique that constructs a polynomial from data points (xi, yi) where the xi need NOT be equally spaced. 

A divided difference table is built recursively:
- First divided difference: [xi, xi+1] = (f(xi+1) – f(xi)) / (xi+1 – xi).
- Higher-order differences use previous divided differences to form an upper triangular table. 

The interpolating polynomial in one line is:
P(x) = f[x0] + f[x0, x1](x – x0) + f[x0, x1, x2](x – x0)(x – x1) + … + f[x0, x1, …, xn](x – x0)(x – x1)…(x – xn-1).
Algorithm: Newton’s Divided Difference

1. Input:
   - Read data points (xi, yi) for i = 0, 1, …, n. 

2. Initialize table:
   - Set dd[i][0] = yi for all i. 

3. Build divided difference table:
   - For k = 1 to n:
       For i = 0 to n - k:
         dd[i][k] = (dd[i+1][k-1] - dd[i][k-1]) / (x_{i+k} - x_i). 

4. Form polynomial:
   - Coefficients are: dd[0][0], dd[0][1], …, dd[0][n].
   - Polynomial:
     P(x) = dd[0][0]
          + dd[0][1](x - x0)
          + dd[0][2](x - x0)(x - x1)
          + … + dd[0][n](x - x0)…(x - x_{n-1})

#### Newton Divided Difference Code
```cpp
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

int main()
{
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    if (!fin)
    {
        cout << "Error opening input file!" << endl;
        return 1;
    }

    int n;
    fin >> n;                 
    int m = n + 1;            

    double x[m];
    double dd[m][m];


    for (int i = 0; i < m; i++)
        fin >> x[i];

    
    for (int i = 0; i < m; i++)
        fin >> dd[i][0];

    
    for (int j = 1; j < m; j++)
    {
        for (int i = 0; i < m - j; i++)
        {
            dd[i][j] =
                (dd[i + 1][j - 1] - dd[i][j - 1]) /
                (x[i + j] - x[i]);
        }
    }

    
    fout << fixed << setprecision(6);
    fout << "Divided Difference Table:\n";

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < m - i; j++)
            fout << setw(12) << dd[i][j] << " ";
        fout << "\n";
    }

    double xp;
    fin >> xp;

    
    double P = dd[0][0];
    double term = 1.0;

    for (int i = 1; i < n; i++)
    {
        term *= (xp - x[i - 1]);
        P += term * dd[0][i];
    }


    double errorFactor = 1.0;
    for (int i = 0; i < m; i++)
        errorFactor *= (xp - x[i]);

    double error = errorFactor * dd[0][n];

    fout << "\nInterpolated value at x = " << xp << " : " << P << "\n";
    fout << "Estimated error : " << error << "\n";

    fin.close();
    fout.close();

    return 0;
}
````
#### Newton Divided Difference Input
```
3
2 3 4 5
0.6931 1.0986 1.3863 1.6094
3.5

````

#### Newton Divided Difference Output
````
Divided Difference Table:
    0.693100     0.405500    -0.058900     0.008867 
    1.098600     0.287700    -0.032300 
    1.386300     0.223100 
    1.609400 

Interpolated value at x = 3.500000 : 1.257175
Estimated error : 0.004987
````
[Back to Table of Contents](#table-of-contents)

---
## Numerical Differentiation
### Differentiation Forward Interpolation Method
#### Differentiation Forward Theory
The Forward Interpolation Differentiation Method is a numerical technique used to approximate the derivative of a function when the function values are known at equally spaced points. The method is based on Newton’s Forward Interpolation formula and is suitable when the value of the derivative is required near the beginning of the data table.

Mathematical / Formula Representation

Newton’s Forward Interpolation formula is:

y = y₀ + uΔy₀ + u(u−1)/2! Δ²y₀ + u(u−1)(u−2)/3! Δ³y₀ + …

where
u = (x − x₀)/h
h = interval between successive x values

Differentiating with respect to x:

dy/dx = (1/h) [ Δy₀ − (1/2)Δ²y₀ + (1/3)Δ³y₀ − (1/4)Δ⁴y₀ + … ]

This formula gives the approximate first derivative at x = x₀.

Algorithm (Short Steps)
1. Read the equally spaced data points (x, y).
2. Construct the forward difference table.
3. Find the step size h.
4. Use the forward differentiation formula.
5. Substitute the forward differences in the formula.
6. Compute the approximate derivative value.

#### Differentiation Forward Code
```cpp
#include <bits/stdc++.h>
using namespace std;

int main()
{
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    int T;
    fin >> T;

    for(int i = 1; i <= T; i++)
    {
    int n;
    fin >> n;

    vector<double> x(n), y(n);

    for (int i = 0; i < n; i++) fin >> x[i];
    for (int i = 0; i < n; i++) fin >> y[i];

    double xp;
    fin >> xp;

    double h = x[1] - x[0]; 
    double p = (xp - x[0]) / h;

    vector<vector<double>> diff(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) diff[i][0] = y[i];

    for (int j = 1; j < n; j++)
        for (int i = 0; i < n - j; i++) diff[i][j] = diff[i + 1][j - 1] - diff[i][j - 1];

    double dydx = 0.0;
    if (n >= 2) dydx += diff[0][1];
    if (n >= 3) dydx += ((2 * p - 1) / 2.0) * diff[0][2];
    if (n >= 4) dydx += ((3 * p * p - 6 * p + 2) / 6.0) * diff[0][3];
    if (n >= 5) dydx += ((4 * p * p * p - 18 * p * p + 22 * p - 6) / 24.0) * diff[0][4];
    dydx /= h;

    fout<<"Test Case " << i << endl;
    fout << "Number of points: " << n << endl;

    fout << fixed << setprecision(3);
    fout<<"x values: ";
    for (int i = 0; i < n; i++) fout << x[i] << " ";
    fout<<endl;
    fout<<"y values: ";
    for (int i = 0; i < n; i++) fout << y[i] << " ";
    fout << "\nStep size(h): " << h << endl;
    fout << "Differentiation point: " << xp << endl;
    fout << "Derivative: " << dydx << endl;

    fout << "\nForward Difference Table (Matrix Form):\n";
    for (int i = 0; i < n; i++)
    {
        fout << setw(4) << x[i] << "\t";
        for (int j = 0; j < n - i; j++)
        {
            double v = diff[i][j];
            if(fabs(v) < 1e-9) v = 0.0;
            fout << setw(12) << v << "\t";
        }
        fout << endl;
    }

    if(i < T) fout << endl;
    }

    fin.close();
    fout.close();

    return 0;
}


````
#### Differentiation Forward Input
```
2
6
1.5 2.0 2.5 3.0 3.5 4.0
3.375 7.0 13.625 24.0 38.875 59.0
1.5

6
1.0 1.2 1.4 1.6 1.8 2.0
0 0.1280 0.5440 1.2960 2.4320 4.00
1.1
````
#### Differentiation Forward Output
```
Test Case 1
Number of points: 6
x values: 1.500 2.000 2.500 3.000 3.500 4.000 
y values: 3.375 7.000 13.625 24.000 38.875 59.000 
Step size(h): 0.500
Differentiation point: 1.500
Derivative: 4.750

Forward Difference Table (Matrix Form):
1.500	       3.375	       3.625	       3.000	       0.750	       0.000	       0.000	
2.000	       7.000	       6.625	       3.750	       0.750	       0.000	
2.500	      13.625	      10.375	       4.500	       0.750	
3.000	      24.000	      14.875	       5.250	
3.500	      38.875	      20.125	
4.000	      59.000	

Test Case 2
Number of points: 6
x values: 1.000 1.200 1.400 1.600 1.800 2.000 
y values: 0.000 0.128 0.544 1.296 2.432 4.000 
Step size(h): 0.200
Differentiation point: 1.100
Derivative: 0.630

Forward Difference Table (Matrix Form):
1.000	       0.000	       0.128	       0.288	       0.048	       0.000	       0.000	
1.200	       0.128	       0.416	       0.336	       0.048	       0.000	
1.400	       0.544	       0.752	       0.384	       0.048	
1.600	       1.296	       1.136	       0.432	
1.800	       2.432	       1.568	
2.000	       4.000	

````
[Back to Table of Contents](#table-of-contents)

---

### Differentiation Backward Interpolation Method
#### Differentiation Backward Theory
The Backward Interpolation Differentiation Method is a numerical technique used to approximate the derivative of a function when the function values are known at equally spaced points. This method is based on Newton’s Backward Interpolation formula and is suitable when the value of the derivative is required near the end of the data table.

Mathematical / Formula Representation

Newton’s Backward Interpolation formula is:

y = yₙ + u∇yₙ + u(u+1)/2! ∇²yₙ + u(u+1)(u+2)/3! ∇³yₙ + …

where
u = (x − xₙ)/h
h = interval between successive x values
∇ = backward difference operator

Differentiating with respect to x:

dy/dx = (1/h) [ ∇yₙ + (1/2)∇²yₙ + (1/3)∇³yₙ + (1/4)∇⁴yₙ + … ]

This formula gives the approximate first derivative at x = xₙ.

Algorithm (Short Steps)
1. Read the equally spaced data points (x, y).
2. Construct the backward difference table.
3. Find the step size h.
4. Use the backward differentiation formula.
5. Substitute the backward differences in the formula.

#### Differentiation Backward Code
```cpp
#include <bits/stdc++.h>
using namespace std;

int main()
{
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    int T;
    fin >> T;

    for(int i = 1; i <= T; i++)
    {
    int n;
    fin >> n;

    vector<double> x(n), y(n);

    for(int i = 0; i < n; i++) fin >> x[i];
    for(int i = 0; i < n; i++) fin >> y[i];

    double xp;
    fin >> xp;

    double h = x[1] - x[0];
    double p = (xp - x[n-1]) / h;

    vector<vector<double>> diff(n, vector<double>(n));
    for(int i = 0; i < n; i++) diff[i][0] = y[i];

    for(int j = 1; j < n; j++)
        for(int i = n-1; i >= j; i--)
            diff[i][j] = diff[i][j-1] - diff[i-1][j-1];

    double dydx = 0.0;

    if(n >= 2) dydx += diff[n-1][1];
    if(n >= 3) dydx += (2*p + 1) * diff[n-1][2] / 2.0;
    if(n >= 4) dydx += (3*p*p + 6*p + 2) * diff[n-1][3] / 6.0;
    if(n >= 5) dydx += (4*p*p*p + 18*p*p + 22*p + 6) * diff[n-1][4] / 24.0;

    dydx /= h;

    fout << "Test Case " << i << endl;
    fout << "Number of points: " << n << endl;

    fout << fixed << setprecision(4);
    fout<<"x values: ";
    for(int i = 0; i < n; i++) fout << x[i] << " ";
    fout << endl;
    fout << "y values: ";
    for(int i = 0; i < n; i++) fout << y[i] << " ";
    fout << "\nStep size(h): " << h << endl;
    fout << "Differentiation point: " << xp << endl;
    fout << "Derivative: " << dydx << endl;

    fout << "\nBackward Difference Table:\n";
    for(int i = 0; i < n; i++)
    {
        fout << setw(4) << x[i] << "\t";
        for(int j = 0; j <= i; j++)
            fout << setw(12) << diff[i][j];
        fout << endl;
    }

    if(i < T) fout << endl;
    }

    return 0;
}

````
#### Differentiation Backward Input
```
2
7
1.0 1.2 1.4 1.6 1.8 2.0 2.2
2.7183 3.3201 4.0552 4.9530 6.0496 7.3891 9.0250
2.2

6
1.0 1.2 1.4 1.6 1.8 2.0
0 0.1280 0.5440 1.2960 2.4320 4.00
1.1
````
#### Differentiation Backward Output
```
Test Case 1
Number of points: 7
x values: 1.0000 1.2000 1.4000 1.6000 1.8000 2.0000 2.2000 
y values: 2.7183 3.3201 4.0552 4.9530 6.0496 7.3891 9.0250 
Step size(h): 0.2000
Differentiation point: 2.2000
Derivative: 9.0214

Backward Difference Table:
1.0000	      2.7183
1.2000	      3.3201      0.6018
1.4000	      4.0552      0.7351      0.1333
1.6000	      4.9530      0.8978      0.1627      0.0294
1.8000	      6.0496      1.0966      0.1988      0.0361      0.0067
2.0000	      7.3891      1.3395      0.2429      0.0441      0.0080      0.0013
2.2000	      9.0250      1.6359      0.2964      0.0535      0.0094      0.0014      0.0001

Test Case 2
Number of points: 6
x values: 1.0000 1.2000 1.4000 1.6000 1.8000 2.0000 
y values: 0.0000 0.1280 0.5440 1.2960 2.4320 4.0000 
Step size(h): 0.2000
Differentiation point: 1.1000
Derivative: 0.6300

Backward Difference Table:
1.0000	      0.0000
1.2000	      0.1280      0.1280
1.4000	      0.5440      0.4160      0.2880
1.6000	      1.2960      0.7520      0.3360      0.0480
1.8000	      2.4320      1.1360      0.3840      0.0480      0.0000
2.0000	      4.0000      1.5680      0.4320      0.0480      0.0000      0.0000
````
[Back to Table of Contents](#table-of-contents)

---
## Curve Fitting: Regression
### Linear Regression Method
#### Linear Regression Theory
The Curve Fitting Linear Regression Method is a numerical technique used to determine a straight-line relationship between two variables based on experimental or observed data. The method fits a linear equation to the data such that the sum of the squares of the deviations between the observed values and the computed values is minimized.

Matrix Representation

The linear curve is assumed in the form:

y = a + b x

For n observed data points, the normal equations are:

Σy = n a + b Σx
Σxy = a Σx + b Σx²

In matrix form:

[ n Σx ] [ a ] = [ Σy ]
[ Σx Σx² ] [ b ] [ Σxy ]

Solving this system gives the coefficients a and b.

Algorithm (Short Steps)
1. Read the given data points (x, y).
2. Assume the linear model y = a + b x.
3. Compute Σx, Σy, Σx², and Σxy.
4. Form the normal equations in matrix form.
5. Solve the equations to obtain a and b.
6. Write the fitted curve using the obtained values.

#### Linear Regression Code
```cpp
#include <bits/stdc++.h>
using namespace std;

int main()
{
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    if(!fin || !fout)
    {
        cout<<"Error"<<endl;
        return 0;
    }

    int T;
    fin>>T;

    for(int i = 1; i <= T; i++)
    {
    int n;
    fin>>n;    
    vector<double> x(n), y(n);

    for(int i = 0; i < n; i++) fin>>x[i];

    for(int i = 0; i < n; i++) fin>>y[i];

    double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0;

    for (int i = 0; i < n; i++) {
        sumX += x[i];
        sumY += y[i];
        sumXY += x[i] * y[i];
        sumX2 += x[i] * x[i];
    }

    double b = (n * sumXY - sumX * sumY) / (n * sumX2 - sumX * sumX);
    double a = (sumY - b * sumX) / n;

    fout << fixed << setprecision(3);

    fout<<"Test Case "<<i<<endl;

    fout<<"Number of points: "<<n;

    fout<<"\nx values: ";
    for(int i=0;i<n;i++) fout << x[i]<<" ";

    fout<<"\ny values: ";
    for(int i=0;i<n;i++) fout << y[i]<<" ";

    fout<<"\nIntercept(a): "<<a;
    fout<<"\nSlope(b): "<<b;
    fout << "\nLinear Regression Equation: y = " << a << " + " << b << "x" << endl;

    if(i < T) fout<<endl;
    }

    fin.close();
    fout.close();

    return 0;
}
```
#### Linear Regression Input
```
3
7
1 2 3 4 5 6 7
3 4 4 5 8 9 10

4
3 9 5 3
8 6 4 2

4
4 7 3 1
6 5 8 3
````
#### Linear Regression Output
```
Test Case 1
Number of points: 7
x values: 1.000 2.000 3.000 4.000 5.000 6.000 7.000 
y values: 3.000 4.000 4.000 5.000 8.000 9.000 10.000 
Intercept(a): 1.143
Slope(b): 1.250
Linear Regression Equation: y = 1.143 + 1.250x

Test Case 2
Number of points: 4
x values: 3.000 9.000 5.000 3.000 
y values: 8.000 6.000 4.000 2.000 
Intercept(a): 4.16l
Slope(b): 0.167
Linear Regression Equation: y = 4.167 + 0.167x

Test Case 3
Number of points: 4
x values: 4.000 7.000 3.000 1.000 
y values: 6.000 5.000 8.000 3.000 
Intercept(a): 4.800
Slope(b): 0.187
Linear Regression Equation: y = 4.800 + 0.187x

````
[Back to Table of Contents](#table-of-contents)

---
### Polynomial Regression Method
#### Polynomial Regression Theory
The Curve Fitting Polynomial Regression Method is a numerical technique used to find a polynomial equation that best fits a given set of experimental or observed data. The method determines the coefficients of the polynomial such that the sum of the squares of the deviations between observed values and calculated values is minimized.

Matrix Representation

The polynomial curve of degree m is assumed as:

y = a₀ + a₁x + a₂x² + … + aₘxᵐ

For n observed data points, the normal equations can be written in matrix form as:

[ n Σx Σx² … Σxᵐ ] [ a₀ ] [ Σy ]
[ Σx Σx² Σx³ … Σxᵐ⁺¹ ] [ a₁ ] [ Σxy ]
[ Σx² Σx³ Σx⁴ … Σxᵐ⁺² ] [ a₂ ] = [ Σx²y ]
[ . . . … . ] [ . ] [ . ]
[ Σxᵐ Σxᵐ⁺¹ Σxᵐ⁺² … Σx²ᵐ ] [ aₘ ] [ Σxᵐy ]

Solving this system gives the polynomial coefficients.

Algorithm (Short Steps)
1. Read the given data points (x, y).
2. Assume a polynomial of degree m.
3. Compute required summations: Σx, Σx², …, Σx²ᵐ and Σy, Σxy, …, Σxᵐy.
4. Form the normal equations in matrix form.
5. Solve the system to obtain the coefficients a₀, a₁, …, aₘ.
6. Write the fitted polynomial curve.

#### Polynomial Regression Code
```cpp
#include <bits/stdc++.h>
using namespace std;

vector<double> gaussianElimination(vector<vector<double>> A, vector<double> b)
{
    int n = A.size();

    for (int i = 0; i < n; i++)
    {
        int pivot = i;
        for (int j = i + 1; j < n; j++)
            if (fabs(A[j][i]) > fabs(A[pivot][i])) pivot = j;

        swap(A[i], A[pivot]);
        swap(b[i], b[pivot]);

        double div = A[i][i];
        for (int j = i; j < n; j++) A[i][j] /= div;
        b[i] /= div;

        for (int j = 0; j < n; j++)
        {
            if (j != i)
            {
                double mul = A[j][i];
                for (int k = i; k < n; k++) A[j][k] -= mul * A[i][k];
                b[j] -= mul * b[i];
            }
        }
    }
    return b;
}

int main()
{
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    int T;
    fin>>T;

    for(int i = 1; i <= T; i++)
    {
    int m, deg;
    fin >> m >> deg;
    
    int n = deg + 1;  

    vector<double> x(m), y(m);

    for (int i = 0; i < m; i++) fin >> x[i];
    for (int i = 0; i < m; i++) fin >> y[i];  

    vector<vector<double>> A(n, vector<double>(n, 0));
    vector<double> B(n, 0);

    vector<double> S(2 * deg + 1, 0);
    for (int i = 0; i < m; i++)
        for (int p = 0; p <= 2 * deg; p++) S[p] += pow(x[i], p);

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) A[i][j] = S[i + j];

    for (int i = 0; i < n; i++)
        for (int k = 0; k < m; k++) B[i] += y[k] * pow(x[k], i);

    vector<double> coeff = gaussianElimination(A, B);

    fout << "Test Case " << i << endl;
    fout << "Data Points: " << m << endl;
    fout << "Polynomial Degree: " << deg << endl;

    fout << "x values: ";
    for (int i = 0; i < m; i++) fout << x[i] << " ";
    fout << "\ny values: ";
    for (int i = 0; i < m; i++) fout << y[i] << " "; 
    fout << endl;

    fout << fixed << setprecision(2);
    fout << "\nCoefficients:\n";
    for (int i = 0; i < n; i++) fout << "a[" << i << "] = " << coeff[i] << endl;

    fout << "\nFitted Polynomial:\n";
    fout << "y = ";
    for (int i = 0; i < n; i++)
    {
    if (i > 0) fout << " + ";
    fout << coeff[i];
    if (i == 1) fout << "x";
    else if (i > 1) fout << "x^" << i;
    }
    fout << endl;

    if(i < T) fout<<endl;
    }

    fin.close();
    fout.close();

    return 0;
}

````
#### Polynomial Regression Input
```
2
5 2
1 2 3 4 5
6 11 18 27 38

5 2
3 4 5 6 7
2.5 3.2 3.8 6.5 11.5
```
#### Polynomial Regression Output
```
Test Case 1
Data Points: 5
Polynomial Degree: 2
x values: 1 2 3 4 5 
y values: 6 11 18 27 38 

Coefficients:
a[0] = 3.00
a[1] = 2.00
a[2] = 1.00

Fitted Polynomial:
y = 3.00 + 2.00x + 1.00x^2

Test Case 2
Data Points: 5
Polynomial Degree: 2
x values: 3.00 4.00 5.00 6.00 7.00 
y values: 2.50 3.20 3.80 6.50 11.50 

Coefficients:
a[0] = 12.43
a[1] = -5.51
a[2] = 0.76

Fitted Polynomial:
y = 12.43 + -5.51x + 0.76x^2

````
[Back to Table of Contents](#table-of-contents)

---
  
### Transcendental Regression Method
#### Transcendental Regression Theory
The Curve Fitting Transcendental Regression Method is used when the relationship between variables is non-polynomial and follows a transcendental form such as exponential or power functions. The method transforms the given nonlinear equation into a linear form using logarithms and then applies linear regression to determine the constants.

Matrix / Linearized Representation

Common transcendental models are:

Exponential model
y = a e^(b x)

Taking logarithm:
ln y = ln a + b x

Let
Y = ln y, A = ln a

Then the linear form becomes:
Y = A + b x

Power model
y = a x^b

Taking logarithm:
ln y = ln a + b ln x

Let
Y = ln y, X = ln x, A = ln a

Then the linear form becomes:
Y = A + b X

After linearization, normal equations of linear regression are applied.

Algorithm (Short Steps)
1. Read the given data points (x, y).
2. Select a suitable transcendental model (exponential or power).
3. Convert the model into linear form using logarithms.
4. Compute required summations for linear regression.
5. Solve the normal equations to obtain constants A and b.
6. Find a by taking antilog of A.
7. Write the fitted transcendental curve.

#### Transcendental Regression Code
```cpp
#include<bits/stdc++.h>
using namespace std;

pair<double, double>PolyRegression(vector<double>&x, vector<double>&y)
{
    int n = x.size();
    double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0;

    for (int i = 0; i < n; i++) {
        sumX += x[i];
        sumY += y[i];
        sumXY += x[i] * y[i];
        sumX2 += x[i] * x[i];
    }

    double B = (n * sumXY - sumX * sumY) / (n * sumX2 - sumX * sumX);
    double A = (sumY - B * sumX) / n;

    return {A, B};
}

int main()
{
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    int T;
    fin >> T;

    for(int i = 1; i <= T; i++)
    {
    int Type;
    fin>>Type;  

    int n;
    fin>>n;    
    vector<double> x(n), y(n);

    for(int i = 0; i < n; i++) fin>>x[i];
    for(int i = 0; i < n; i++) fin>>y[i];

    fout << fixed << setprecision(3);
    fout << "Test Case " << i << endl;
    fout << "Equation Type: " << Type << endl;
    fout << "Number of points: " << n << endl;

    fout<<"x values: ";
    for(int i=0;i<n;i++) fout << x[i]<<" ";

    fout<<endl;

    fout<<"y values: ";
    for(int i=0;i<n;i++) fout << y[i]<<" ";

    vector<double>X(n), Y(n);

    pair<double,double> res;
    double A, B, a, b, k;

    switch(Type)
    {
    case 1:
    {
        for(int i = 0; i < n; i++)
        {
            Y[i] = log(y[i]);
            X[i] = log(x[i]);
        }

        res = PolyRegression(X, Y);
        A = res.first;
        B = res.second;

        a = exp(A);
        b = B;

        fout << endl;
        fout << "a: " << a << endl;
        fout << "b: " << b << endl;

        fout << "Transcendental Equation: ";
        fout << "y = " << a << "*x^" << "(" << b << ")" << endl;
        break;
    }

    case 2:
    {
        for(int i = 0; i < n; i++)
        {
            Y[i] = log(y[i]);
            X[i] = x[i];
        }

        res = PolyRegression(X, Y);
        A = res.first;
        B = res.second;

        a = exp(A);
        k = B;

        fout << endl;
        fout << "a: " << a << endl;
        fout << "k: " << k << endl;

        fout << "Transcendental Equation: ";
        fout << "y = " << a << "*e^(" << k << "*x)" << endl;
        break;
    }

    case 3:
    {
        for(int i = 0; i < n; i++)
        {
            Y[i] = y[i];
            X[i] = exp(x[i] / 4.0);
        }

        res = PolyRegression(X, Y);
        a = res.first;
        b = res.second;

        fout << endl;
        fout << "a: " << a << endl;
        fout << "b: " << b << endl;

        fout << "Transcendental Equation: ";
        fout << "y = " << a << " + " << b << "*e^(x/4.0)" << endl;
        break;
    }

    default:
    {
        fout << "Invalid Type." << endl;
        break;
    }
    }

    if(i < T) fout<<endl;
    }

    fin.close();
    fout.close();

    return 0;
}

````
#### Transcendental Regression Input
```
3
1
4
2 3 4 5
8 27 64 125

2
5
1 2 3 4 5
2.7 7.4 20.1 54.6 148.4

3
5
1 2 3 4 5
50 80 96 120 145
```
#### Transcendental Regression Output
```
Test Case 1
Equation Type: 1
Number of points: 4
x values: 2.000 3.000 4.000 5.000 
y values: 8.000 27.000 64.000 125.000 
a: 1.000
b: 3.000
Transcendental Equation: y = 1.000*x^(3.000)

Test Case 2
Equation Type: 2
Number of points: 5
x values: 1.000 2.000 3.000 4.000 5.000 
y values: 2.700 7.400 20.100 54.600 148.400 
a: 0.996
k: 1.001
Transcendental Equation: y = 0.996*e^(1.001*x)

Test Case 3
Equation Type: 3
Number of points: 5
x values: 1.000 2.000 3.000 4.000 5.000 
y values: 50.000 80.000 96.000 120.000 145.000 
a: 5.749
b: 41.059
Transcendental Equation: y = 5.749 + 41.059*e^(x/4.0)
```
[Back to Table of Contents](#table-of-contents)
