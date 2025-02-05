#include <iostream>
#include "math/Vector2D.hpp"

using namespace std;

void printCoords(Vector2D& v) { 
	cout << "(x,y) = (" << v.getX() << "," << v.getY() << ")\n";
}

int main() {
	Vector2D v1;
	Vector2D v2(3.5, -1.0);

	printCoords(v1); printCoords(v2);
}
