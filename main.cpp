// PointInTriangle.cpp

#include <math.h>
#include <cstdint>
#include <stdlib.h>
#include <iostream>
#include <chrono>

#define EPSYLON 0.000001
#define NUMBER_OF_TESTS 10000000

struct Point {
	double X;
	double Y;

	Point() {};

	Point(double x, double y) {
		X = x;
		Y = y;
	}
};

Point points_tested[NUMBER_OF_TESTS];
Point points_a[NUMBER_OF_TESTS];
Point points_b[NUMBER_OF_TESTS];
Point points_c[NUMBER_OF_TESTS];
bool barycentric_results[NUMBER_OF_TESTS];
bool orientation_results[NUMBER_OF_TESTS];

//The first approach
bool PointLiesInTriangle_Barycentric(Point& A, Point& B, Point &C, Point& P) {
	double BY_CY = B.Y - C.Y;
	double CX_BX = C.X - B.X;
	double BXCY_CXBY = B.X * C.Y - C.X * B.Y;
	double CY_AY = C.Y - A.Y;
	double AX_CX = A.X - C.X;
	double CXAY_AXCY = C.X * A.Y - A.X * C.Y;
	double AY_BY = A.Y - B.Y;
	double BX_AX = B.X - A.X;
	double AXBY_BXAY = A.X * B.Y - B.X * A.Y;

	// Calculate the determinant of barycentric coordinates matrix
	double D = (A.X * (BY_CY)+B.X * (CY_AY)+C.X * (AY_BY));
	if (fabs(D) < EPSYLON)
		return false;
	D = 1 / D;

	double A1 = (P.X * (BY_CY) + P.Y * (CX_BX) + (BXCY_CXBY)) * D;
	double A2 = (P.X * (CY_AY) + P.Y * (AX_CX) + (CXAY_AXCY)) * D;
	double A3 = (P.X * (AY_BY) + P.Y * (BX_AX) + (AXBY_BXAY)) * D;

	return  A1 >= 0 && A2 >= 0 && A3 >= 0 && A1 + A2 + A3 <= 1 + EPSYLON;
}

enum class Orientation
{
	Clockwise,
	Counter_Clockwise,
	Collinear
};

//Note that this function id only used when generating the tests, not when the results are actually computed!!!
Orientation GetOrientation(double x1, double y1, double x2, double y2, double x3, double y3) {
	double val = (y2 - y1) * (x3 - x2) - (x2 - x1) * (y3 - y2);
	if (val < EPSYLON && val > -EPSYLON)
		return Orientation::Collinear;
	else
		return val < 0 ? Orientation::Counter_Clockwise : Orientation::Clockwise;
}

//The second (and the best) approach
static bool PointLiesInTriangle_Orientation(Point& A, Point& B, Point& C, Point& P) {
	if (fabs(A.X - P.X) < EPSYLON && fabs(A.Y - P.Y) < EPSYLON ||
		fabs(B.X - P.X) < EPSYLON && fabs(B.Y - P.Y) < EPSYLON ||
		fabs(C.X - P.X) < EPSYLON && fabs(C.Y - P.Y) < EPSYLON)
		return true;

	double val1 = (A.Y - P.Y) * (B.X - A.X) - (A.X - P.X) * (B.Y - A.Y);
	if (fabs(val1) < EPSYLON) {
		//if (val1 < EPSYLON && val1 > -EPSYLON) {
		return ((P.X >= A.X && P.X <= B.X) || (P.X >= B.X && P.X <= A.X)) &&
			   ((P.Y >= A.Y && P.Y <= B.Y) || (P.Y >= B.Y && P.Y <= A.Y));
	}
	Orientation orientation_first = val1 < 0 ? Orientation::Counter_Clockwise : Orientation::Clockwise;

	double val2 = (B.Y - P.Y) * (C.X - B.X) - (B.X - P.X) * (C.Y - B.Y);
	if (fabs(val2) < EPSYLON) {
		//if(val2 < EPSYLON && val2 > -EPSYLON) {
		return ((P.X >= B.X && P.X <= C.X) || (P.X >= C.X && P.X <= B.X)) &&
			   ((P.Y >= B.Y && P.Y <= C.Y) || (P.Y >= C.Y && P.Y <= B.Y));
	}
	Orientation orientation_central = val2 < 0 ? Orientation::Counter_Clockwise : Orientation::Clockwise;

	double val3 = (C.Y - P.Y) * (A.X - C.X) - (C.X - P.X) * (A.Y - C.Y);
	if (fabs(val3) < EPSYLON) {
		//if (val3 < EPSYLON && val3 > -EPSYLON) {
		return ((P.X >= C.X && P.X <= A.X) || (P.X >= A.X && P.X <= C.X)) &&
			   ((P.Y >= C.Y && P.Y <= A.Y) || (P.Y >= A.Y && P.Y <= C.Y));
	}
	Orientation orientation_third = val3 < 0 ? Orientation::Counter_Clockwise : Orientation::Clockwise;

	return orientation_first == orientation_central &&
		   orientation_central == orientation_third;
}

int main() {
	std::chrono::time_point<std::chrono::steady_clock> milliseconds_for_random = std::chrono::steady_clock::now();
	srand(std::chrono::duration_cast<std::chrono::milliseconds>(milliseconds_for_random.time_since_epoch()).count() % (0x100000000));
	bool overall_success = true;

	for (int32_t i = 0; i < NUMBER_OF_TESTS; i++) {
		int32_t x1 = rand() % 1000;
		int32_t x2 = rand() % 1000;
		int32_t x3 = rand() % 1000;
		int32_t y1 = rand() % 1000;
		int32_t y2 = rand() % 1000;
		int32_t y3 = rand() % 1000;
		while (GetOrientation(x1, y1, x2, y2, x3, y3) == Orientation::Collinear) {
			x1 = rand() % 1000;
			x2 = rand() % 1000;
			x3 = rand() % 1000;
			y1 = rand() % 1000;
			y2 = rand() % 1000;
			y3 = rand() % 1000;
		}
		int32_t tested_x = rand() % 1000;
		int32_t tested_y = rand() % 1000;
		points_a[i] = Point(x1, y1);
		points_b[i] = Point(x2, y2);
		points_c[i] = Point(x3, y3);
		points_tested[i] = Point(tested_x, tested_y);
	}

	std::chrono::time_point<std::chrono::steady_clock> milliseconds_start = std::chrono::steady_clock::now();
	for (int32_t i = 0; i < NUMBER_OF_TESTS; i++)
		barycentric_results[i] = PointLiesInTriangle_Barycentric(points_a[i], points_b[i], points_c[i], points_tested[i]);
	auto delta1 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - milliseconds_start);

	milliseconds_start = std::chrono::steady_clock::now();
	for (int32_t i = 0; i < NUMBER_OF_TESTS; i++)
		orientation_results[i] = PointLiesInTriangle_Orientation(points_a[i], points_b[i], points_c[i], points_tested[i]);
	auto delta2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - milliseconds_start);

	for (int32_t i = 0; i < NUMBER_OF_TESTS; i++)
		if (barycentric_results[i] != orientation_results[i]) {
			overall_success = false;
			break;
		}

	if (overall_success)
		std::cout << "OVERALL : SUCCESS" << std::endl;
	else
		std::cout << "OVERALL : FAILURE" << std::endl;
	std::cout << "barycentic coordinates time " << static_cast<double>(delta1.count()) << std::endl;
	std::cout << " orientation time " << static_cast<double>(delta2.count()) << std::endl;
	return 0;
}
