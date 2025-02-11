#include <iostream>
#include "math/Vector2D.hpp"

using namespace std;

void printCoords(const Vector2D& v) {
    cout << "Vector: " << v << endl;
}

void testVector2D() {
    // Default Constructor
    Vector2D v1;
    Vector2D v2(3.5, -1.0);

    cout << "Testing Constructors:\n";
    printCoords(v1);
    printCoords(v2);

    // Arithmetic Operators
    cout << "\nTesting Arithmetic Operators:\n";
    Vector2D v3 = v1 + v2;
    cout << "v1 + v2 = " << v3 << endl;

    Vector2D v4 = v2 - v1;
    cout << "v2 - v1 = " << v4 << endl;

    Vector2D v5 = v2 * 2.0;
    cout << "v2 * 2 = " << v5 << endl;

    Vector2D v6 = v2 / 2.0;
    cout << "v2 / 2 = " << v6 << endl;

    // Compound Assignment Operators
    cout << "\nTesting Compound Assignment Operators:\n";
    Vector2D v7 = v2;
    v7 += v1;
    cout << "v7 += v1 → " << v7 << endl;

    v7 -= v2;
    cout << "v7 -= v2 → " << v7 << endl;

    v7 *= 3;
    cout << "v7 *= 3 → " << v7 << endl;

    v7 /= 2;
    cout << "v7 /= 2 → " << v7 << endl;

    // Comparison Operators
    cout << "\nTesting Comparison Operators:\n";
    cout << "v1 == v2 → " << (v1 == v2 ? "true" : "false") << endl;
    cout << "v1 != v2 → " << (v1 != v2 ? "true" : "false") << endl;
}

int main() {
    testVector2D();
    return 0;
}

