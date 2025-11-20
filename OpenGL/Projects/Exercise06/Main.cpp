#include <glad/glad.h>
#include <GLFW/glfw3.h>

#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>

// IMPORTANT: If you are using Mac, kindly change this to the full path of the test.txt file
#define INPUT_FILE "test.txt"

/**
 * Point struct
 */
struct Point {
	/**
	 * X-coordinate of the point
	 */
	float x;

	/**
	 * Y-coordinate of the point
	 */
	float y;

	Point(float newX, float newY)
	{
		x = newX;
		y = newY;
	}

	Point()
		: Point(0.0f, 0.0f)
	{
	}

	Point operator-(const Point& other) const {
        return Point(x - other.x, y - other.y);
    }

	float cross(const Point& b) const {
        return x * b.y - y * b.x;
    }
};

/**
 * Triangle struct
 */
struct Triangle {
	/**
	 * First point
	 */
	Point a;

	/**
	 * Second point
	 */
	Point b;

	/**
	 * Third point
	 */
	Point c;

	Triangle(Point newA, Point newB, Point newC)
	: a(newA), b(newB), c(newC)
    {}
};

/**
 * Polygon struct
 */
struct Polygon {
	/**
	 * Point list
	 */
	std::vector<Point> outline;

	/**
	 * Holes list (each hole is a list of points describing the outline of the hole)
	 */
	std::vector< std::vector<Point> > holes;
};

bool isPointInsideTriangle(const Point& A, const Point& B, const Point& C, const Point& P) {
    float cross1 = (B - A).cross(P - A);
    float cross2 = (C - B).cross(P - B);
    float cross3 = (A - C).cross(P - C);
    
    return (cross1 >= 0 && cross2 >= 0 && cross3 >= 0) ||  (cross1 <= 0 && cross2 <= 0 && cross3 <= 0);
}

/**
 * @brief Splits the specified polygon into a set of triangles
 * @param[in] polygon Polygon 
 */
std::vector<Triangle> Triangulate(const Polygon& polygon) {
	std::vector<Triangle> ret;
	std::vector<Point> workingOutline = polygon.outline;

	while (workingOutline.size() > 3) {
		// earclip
		for (int i = 0; i < workingOutline.size(); i++) {
			int prevIndex = (i - 1 + workingOutline.size()) % workingOutline.size();
			int nextIndex = (i + 1) % workingOutline.size();
			
			const Point& vPrev = workingOutline[prevIndex];
			const Point& vCurrent = workingOutline[i];
			const Point& vNext = workingOutline[nextIndex];

			// check if E is convex (cross product > 0)
			float cross = (vCurrent - vPrev).cross(vNext - vCurrent);

			// check if no other vertices inside triangle
            bool hasPointsInside = false;
            for (int j = 0; j < workingOutline.size(); j++) {
                if (j == i || j == prevIndex || j == nextIndex) continue;
                
                if (isPointInsideTriangle(vPrev, vCurrent, vNext, workingOutline[j])) {
                    hasPointsInside = true;
                    break;
                }
            }

			// if both YES, then clip ear (add these 3 vertices as a triangle to the ret, remove vCurrent from polygon outline)
			if (cross > 0 && !hasPointsInside) {
				ret.push_back(Triangle(vPrev, vCurrent, vNext));
				workingOutline.erase(workingOutline.begin() + i);
				break;
			}
		}
	}

	if (workingOutline.size() == 3) {
        ret.push_back(Triangle(workingOutline[0], workingOutline[1], workingOutline[2]));
    }
	return ret;
}

// ----------------------------------------------------------------------------------
// DO NOT TOUCH ANY CODE BEYOND THIS POINT
// ----------------------------------------------------------------------------------
#include <cstddef>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <random>
#include <string>

const char* vertShaderSource = R"(#version 330
layout(location = 0) in vec3 vertexPosition;
layout(location = 1) in vec3 vertexColor;
out vec3 outVertexColor;
uniform mat4 mvpMatrix;
void main() {
    gl_Position = mvpMatrix * vec4(vertexPosition, 1.0);
	outVertexColor = vertexColor;
})";
const char* fragShaderSource = R"(#version 330
in vec3 outVertexColor;
out vec4 fragColor;
void main() {
	fragColor = vec4(outVertexColor, 1.0);
})";
const int WINDOW_WIDTH = 800;
const int WINDOW_HEIGHT = 600;
struct AppData {
	std::vector<Polygon> testCases;
	int currentTestCase;
	GLuint vbo;
	size_t numVertices;
	bool showTriangulation;
};
struct Vertex {
	GLfloat x, y, z;
	GLfloat r, g, b;
};
struct Matrix4x4 {
	GLfloat values[16];
};

void RefreshScene(AppData* data);

Matrix4x4 operator*(const Matrix4x4& a, const Matrix4x4& b) {
	Matrix4x4 ret;
	ret.values[0] = a.values[0] * b.values[0] + a.values[4] * b.values[1] + a.values[8] * b.values[2] + a.values[12] * b.values[3];
	ret.values[4] = a.values[0] * b.values[4] + a.values[4] * b.values[5] + a.values[8] * b.values[6] + a.values[12] * b.values[7];
	ret.values[8] = a.values[0] * b.values[8] + a.values[4] * b.values[9] + a.values[8] * b.values[10] + a.values[12] * b.values[11];
	ret.values[12] = a.values[0] * b.values[12] + a.values[4] * b.values[13] + a.values[8] * b.values[14] + a.values[12] * b.values[15];
	ret.values[1] = a.values[1] * b.values[0] + a.values[5] * b.values[1] + a.values[9] * b.values[2] + a.values[13] * b.values[3];
	ret.values[5] = a.values[1] * b.values[4] + a.values[5] * b.values[5] + a.values[9] * b.values[6] + a.values[13] * b.values[7];
	ret.values[9] = a.values[1] * b.values[8] + a.values[5] * b.values[9] + a.values[9] * b.values[10] + a.values[13] * b.values[11];
	ret.values[13] = a.values[1] * b.values[12] + a.values[5] * b.values[13] + a.values[9] * b.values[14] + a.values[13] * b.values[15];
	ret.values[2] = a.values[2] * b.values[0] + a.values[6] * b.values[1] + a.values[10] * b.values[2] + a.values[14] * b.values[3];
	ret.values[6] = a.values[2] * b.values[4] + a.values[6] * b.values[5] + a.values[10] * b.values[6] + a.values[14] * b.values[7];
	ret.values[10] = a.values[2] * b.values[8] + a.values[6] * b.values[9] + a.values[10] * b.values[10] + a.values[14] * b.values[11];
	ret.values[14] = a.values[2] * b.values[12] + a.values[6] * b.values[13] + a.values[10] * b.values[14] + a.values[14] * b.values[15];
	ret.values[3] = a.values[3] * b.values[0] + a.values[7] * b.values[1] + a.values[11] * b.values[2] + a.values[15] * b.values[3];
	ret.values[7] = a.values[3] * b.values[4] + a.values[7] * b.values[5] + a.values[11] * b.values[6] + a.values[15] * b.values[7];
	ret.values[11] = a.values[3] * b.values[8] + a.values[7] * b.values[9] + a.values[11] * b.values[10] + a.values[15] * b.values[11];
	ret.values[15] = a.values[3] * b.values[12] + a.values[7] * b.values[13] + a.values[11] * b.values[14] + a.values[15] * b.values[15];

	return ret;
}

GLuint CreateShader(const GLuint& type, const std::string& source) {
	GLuint shader = glCreateShader(type);

	const char* sourceCStr = source.c_str();
	GLint sourceLen = source.size();
	glShaderSource(shader, 1, &sourceCStr, &sourceLen);
	glCompileShader(shader);

	GLint compileStatus;
	glGetShaderiv(shader, GL_COMPILE_STATUS, &compileStatus);
	if (compileStatus == GL_FALSE) {
		char infoLog[512];
		GLsizei infoLogLen = sizeof(infoLog);
		glGetShaderInfoLog(shader, infoLogLen, &infoLogLen, infoLog);

		std::string errorMsg;
		if (type == GL_VERTEX_SHADER) {
			errorMsg += std::string("Failed to compile vertex shader!\n");
		}
		else if (type == GL_FRAGMENT_SHADER) {
			errorMsg += std::string("Failed to compile fragment shader!\n");
		}
		else {
			errorMsg += std::string("Failed to compile shader!\n");
		}
		errorMsg += std::string(infoLog);

		std::cout << errorMsg << std::endl;
	}

	return shader;
}

GLuint CreateShaderProgramFromSource(const std::string& vertexShaderSource, const std::string& fragmentShaderSource) {
	GLuint vsh = CreateShader(GL_VERTEX_SHADER, vertexShaderSource);
	GLuint fsh = CreateShader(GL_FRAGMENT_SHADER, fragmentShaderSource);

	GLuint program = glCreateProgram();
	glAttachShader(program, vsh);
	glAttachShader(program, fsh);
	glLinkProgram(program);

	GLint linkStatus;
	glGetProgramiv(program, GL_LINK_STATUS, &linkStatus);
	if (linkStatus != GL_TRUE) {
		char infoLog[512];
		GLsizei infoLogLen = sizeof(infoLog);
		glGetProgramInfoLog(program, infoLogLen, &infoLogLen, infoLog);
		throw std::runtime_error(std::string("program link error: ") + infoLog);
		return 0;
	}

	glDetachShader(program, vsh);
	glDetachShader(program, fsh);
	glDeleteShader(vsh);
	glDeleteShader(fsh);

	return program;
}

void KeyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
	AppData* appData = reinterpret_cast<AppData*>(glfwGetWindowUserPointer(window));
	if (action == GLFW_PRESS) {
		if (key == GLFW_KEY_LEFT) {
			appData->currentTestCase = ((appData->currentTestCase - 1) + appData->testCases.size()) % appData->testCases.size();
		}
		else if (key == GLFW_KEY_RIGHT) {
			appData->currentTestCase = (appData->currentTestCase + 1) % appData->testCases.size();
		}
		else if (key == GLFW_KEY_SPACE) {
			appData->showTriangulation = !appData->showTriangulation;
		}
		RefreshScene(appData);
	}
}

Matrix4x4 CreateIdentity() {
	Matrix4x4 ret = {};
	ret.values[0] = 1.0f; ret.values[5] = 1.0f; ret.values[10] = 1.0f; ret.values[15] = 1.0f;
	return ret;
}

Matrix4x4 CreateOrtho(float left, float right, float bottom, float top, float near, float far) {
	Matrix4x4 ret =
	{
		2.0f / (right - left), 0.0f, 0.0f, 0.0f,
		0.0f, 2.0f / (top - bottom), 0.0f, 0.0f,
		0.0f, 0.0f, -2.0f / (far - near), 0.0f,
		-(right + left) / (right - left), -(top + bottom) / (top - bottom), -(far + near) / (far - near), 1.0f
	};

	return ret;
}

void AppendLineVertices(const Point& p0, const Point& p1, float width, std::vector<Vertex>& vertices, float zOrder = 0.0f, float r = 1.0f, float g = 1.0f, float b = 1.0f) {
	float vx = p1.x - p0.x, vy = p1.y - p0.y;
	float vm = vx * vx + vy * vy;
	if (fabsf(vm) > 1e-9f) {
		vm = sqrtf(vm);
		vx = vx / vm; vy = vy / vm;
	}
	float nx = -vy, ny = vx;

	Vertex v0 = { p0.x + nx * width / 2.0f, p0.y + ny * width / 2.0f, zOrder, r, g, b };
	Vertex v1 = { p0.x - nx * width / 2.0f, p0.y - ny * width / 2.0f, zOrder, r, g, b };
	Vertex v2 = { p1.x - nx * width / 2.0f, p1.y - ny * width / 2.0f, zOrder, r, g, b };
	Vertex v3 = { p1.x + nx * width / 2.0f, p1.y + ny * width / 2.0f, zOrder, r, g, b };

	vertices.push_back( v0 ); vertices.push_back( v1 ); vertices.push_back( v2 );
	vertices.push_back( v2 ); vertices.push_back( v3 ); vertices.push_back( v0 );
}

void AppendCircleVertices(const Point& center, float radius, std::vector<Vertex>& vertices, float zOrder = 0.0f, float r = 1.0f, float g = 1.0f, float b = 1.0f) {
	int numSections = 180;
	float anglePerSection = 360.0f / numSections * M_PI / 180.0f;
	for (int i = 1; i <= numSections; ++i) {
		float x0 = center.x + radius * cosf(anglePerSection * (i - 1));
		float y0 = center.y + radius * sinf(anglePerSection * (i - 1));
		float x1 = center.x + radius * cosf(anglePerSection * i);
		float y1 = center.y + radius * sinf(anglePerSection * i);
		vertices.push_back( { center.x, center.y, zOrder, r, g, b } );
		vertices.push_back( { x0, y0, zOrder, r, g, b } );
		vertices.push_back( { x1, y1, zOrder, r, g, b } );
	}
}

void RefreshScene(AppData* data) {
	std::vector<Vertex> vertices;

	Polygon& polygon = data->testCases[data->currentTestCase];
	for (size_t i = 0; i < polygon.outline.size(); ++i) {
		AppendLineVertices(polygon.outline[i], polygon.outline[(i + 1) % polygon.outline.size()], 2.0f, vertices, 0.0f, 1.0f, 1.0f, 1.0f);
	}
	for (size_t i = 0; i < polygon.holes.size(); ++i) {
		std::vector<Point>& points = polygon.holes[i];
		for (size_t j = 0; j < points.size(); ++j) {
			AppendLineVertices(points[j], points[(j + 1) % points.size()], 2.0f, vertices, 0.0f, 1.0f, 1.0f, 1.0f);
		}
	}

	if (data->showTriangulation) {
		std::vector<Triangle> triangles = Triangulate(data->testCases[data->currentTestCase]);
		for (size_t i = 0; i < triangles.size(); ++i) {
			AppendLineVertices(triangles[i].a, triangles[i].b, 2.0f, vertices, 0.0f, 1.0f, 0.0f, 0.0f);
			AppendLineVertices(triangles[i].b, triangles[i].c, 2.0f, vertices, 0.0f, 1.0f, 0.0f, 0.0f);
			AppendLineVertices(triangles[i].c, triangles[i].a, 2.0f, vertices, 0.0f, 1.0f, 0.0f, 0.0f);
		}
	}

	glBindBuffer(GL_ARRAY_BUFFER, data->vbo);
	glBufferData(GL_ARRAY_BUFFER, sizeof(Vertex) * vertices.size(), vertices.data(), GL_STATIC_DRAW);
	data->numVertices = vertices.size();
}

int main(int argc, char* argv[]) {
	srand(time(nullptr));

	AppData appData = {};
	appData.showTriangulation = true;
	std::ifstream file(INPUT_FILE);
	if (file.fail())
	{
		std::cout << "Failed to read file test.txt" << std::endl;
		return 1;
	}

	size_t numTestCases;
	file >> numTestCases;
	for (size_t testCaseNumber = 0; testCaseNumber < numTestCases; ++testCaseNumber) {
		Polygon polygon;

		size_t numBoundaries;
		file >> numBoundaries;
		for (size_t i = 0; i < numBoundaries; ++i) {
			size_t numPoints;
			file >> numPoints;

			std::vector<Point> points;
			for (size_t j = 0; j < numPoints; ++j) {
				Point point;
				file >> point.x >> point.y;
				points.push_back(point);
			}

			if (i > 0) {
				// Hole
				polygon.holes.push_back(points);
			}
			else {
				polygon.outline = points;
			}
		}

		appData.testCases.push_back(polygon);
	}

	if (glfwInit() == GLFW_FALSE) {
		std::cerr << "Cannot initialize GLFW!" << std::endl;
		return -1;
	}

	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);
	glfwWindowHint(GLFW_SAMPLES, 4);

	GLFWwindow* window = glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "Triangulation", nullptr, nullptr);
	if (!window) {
		std::cerr << "Cannot create window.";
		return -1;
	}

	glfwMakeContextCurrent(window);
	gladLoadGLLoader((GLADloadproc)glfwGetProcAddress);

	glfwSetKeyCallback(window, KeyCallback);

	GLuint vbo, vao;
	glGenBuffers(1, &vbo);
	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), reinterpret_cast<void*>(offsetof(Vertex, x)));
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), reinterpret_cast<void*>(offsetof(Vertex, r)));
	glBindVertexArray(0);

	GLuint program = CreateShaderProgramFromSource(vertShaderSource, fragShaderSource);

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_MULTISAMPLE);

	appData.vbo = vbo;
	glfwSetWindowUserPointer(window, &appData);

	RefreshScene(&appData);

	Matrix4x4 projMatrix = CreateOrtho(0.0f, WINDOW_WIDTH * 1.0f, 0.0f, WINDOW_HEIGHT * 1.0f, -100.0f, 100.0f);
	float prevTime = glfwGetTime(), speed = 200.0f;
	while (!glfwWindowShouldClose(window)) {
		float currentTime = glfwGetTime();
		float deltaTime = currentTime - prevTime;
		prevTime = currentTime;

		glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glUseProgram(program);

		glBindVertexArray(vao);
		glUniformMatrix4fv(glGetUniformLocation(program, "mvpMatrix"), 1, GL_FALSE, projMatrix.values);
		glDrawArrays(GL_TRIANGLES, 0, appData.numVertices);

		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	glDeleteProgram(program); program = 0;
	glDeleteVertexArrays(1, &vao); vao = 0;
	glDeleteBuffers(1, &vbo); vbo = 0;
	glfwDestroyWindow(window);window = nullptr;
	glfwTerminate();

	return 0;
}
