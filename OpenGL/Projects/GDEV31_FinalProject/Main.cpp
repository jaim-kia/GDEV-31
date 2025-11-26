#include <glad/glad.h>
#include <GLFW/glfw3.h>

#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <cstdlib>

/**
 * @brief Path to the input file
 * 
 * IMPORTANT: If you are using Mac, kindly change this to the full path of the test.txt file
 */
#define INPUT_FILE "test.txt"

/**
 * @brief 2D point representation
 */
struct Point {
    /** X-coordinate */
    float x;

    /** Y-coordinate */
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
 * @brief Represents a Voronoi cell for a given site
 */
struct Cell {
    /** Site position */
    Point site;

    /** Vertices forming the boundary of the cell */
    std::vector<Point> vertices;

    /** Red edge color */
    float edgeColor[3] = {1.0f, 0.0f, 0.0f};

	/** Cell color */
	float cellColor[3] = {0.0f, 0.0f, 0.0f};
};

/**
 * @brief Site group used in test cases
 */
struct Sites {
    std::vector<Point> sitelist;
};

float eucDist (Point& a, Point& b) {
	float xDiff = b.x - a.x;
	float yDiff = b.y - a.y;
	float dist = sqrt((xDiff*xDiff) - (yDiff*yDiff));
	return dist;
}

// takes in a polygon and a bisector
// returns a polygon cut or clipped down by the bisector
// tldr: gets the norm of the bisector in order to determine which side to keep (+ norm)
// technically a cell builder 
std::vector<Point> ClipPolygonByLine(const std::vector<Point>& polygon, Point line_p0, Point line_p1) {
    if (polygon.empty()) {
        return {};
    }
    std::vector<Point> clippedPolygon;

    // Calculate normal vector (perpendicular to the bisector)
	// since this is the perp, it points counter-clockwise i.e. away from the site
    float norm_x = -(line_p1.y - line_p0.y);
    float norm_y = (line_p1.x - line_p0.x);

	// goes through all vertices in the polygon and see which ones to retain
    for (int i = 0; i < polygon.size(); ++i) {

		// get the edge of the polygon
        Point current = polygon[i];
        Point next = polygon[(i + 1) % polygon.size()];

        // Calculate signed distance from line creates a vector from 
		// current point to the p0 of the biesector and dot product to the norm
		// we want to check if the created vector is facing the same direction as the norm
        float current_dist = norm_x * (current.x - line_p0.x) + norm_y * (current.y - line_p0.y);
        float next_dist = norm_x * (next.x - line_p0.x) + norm_y * (next.y - line_p0.y);
        
		// meaning site is facing away from the bisector
		// we want all the points that are facing from the norm hence <= 0
        bool currentInside = current_dist <= 0; 
        bool nextInside = next_dist <= 0;

        if (currentInside) {
			// accept when the current point is inside
            clippedPolygon.push_back(current);
        }

		// in case the other point of the edge is not on the side of the bisector
        if (currentInside != nextInside) {
            // Calculate intersection point of the edge and the bisector
			// because it intersects, its distance from the edge is 0.
			// This will be a new vertex that will be part of the cell
			// let t be the ratio between the edge and where the bisector intersects
			// current + t * (next - current) = 0
			//  t * (next - current) = -current
			//  t = -current / (next - current)
			//  t = -current / -(current - next)
			//  t = current / (current - next)
            float t = current_dist / (current_dist - next_dist);  

            Point intersection;
			// get the x coordinate based on the calculated ratio t on the edge
            intersection.x = current.x + t * (next.x - current.x);
			// get the y coordinate based on the calculated ratio t on the edge
            intersection.y = current.y + t * (next.y - current.y);

			// add the new point as a vertex in the polygon
            clippedPolygon.push_back(intersection);
        }

		// Ignore the points if both are outside the bisector
    }
    
    return clippedPolygon;
}

std::vector<Cell> VoronoiDiagram(std::vector<Point>& sites) {
    std::vector<Cell> voronoiCells;

    float min_x = -100000.0f;
    float max_x = 100000.0f;
    float min_y = -100000.0f;
    float max_y = 100000.0f;

    for (int i = 0; i < sites.size(); ++i) {
        Cell cell;
        cell.site = sites[i];
        
        // Start with the bounding box
        cell.vertices = {
            Point(min_x, min_y),
            Point(min_x, max_y),
            Point(max_x, max_y),
            Point(max_x, min_y)
        };

        // Clip against bisectors with all other sites
        for (int j = 0; j < sites.size(); ++j) {
			// check if the bisectors cut out new cells
            if (i == j) continue;

            Point& p_i = sites[i];
            Point& p_j = sites[j];
            
            // Calculate midpoint
            float mid_x = (p_i.x + p_j.x) / 2.0f;
            float mid_y = (p_i.y + p_j.y) / 2.0f;

            // Calculate direction vector from p_i to p_j
            float dir_x = p_j.x - p_i.x;
            float dir_y = p_j.y - p_i.y;
            
            // The bisector is perpendicular to the direction vector
			// extend them with a bigger number to make it into a line on both sides
			// stretching them based on the direction with the mid point as the anchor point
            Point bisector_p0(
                mid_x - dir_y *1000.0f,  
                mid_y + dir_x * 1000.0f
            );
            Point bisector_p1(
                mid_x + dir_y * 1000.0f,   
                mid_y - dir_x * 1000.0f
            );

            // create cells given the calculated bisector
			// start with the "box" made at the start
            cell.vertices = ClipPolygonByLine(cell.vertices, bisector_p0, bisector_p1);
            
            // // Early exit if polygon becomes empty
            // if (cell.vertices.empty()) break;
        }

		// set color of the cell based on the order; make them unique with randomness
		float random_val = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
		float random_color_val = 0 + random_val * (1 - 0);
		cell.cellColor[0] += random_color_val;

		random_val = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
		random_color_val = 0 + random_val * (1 - 0);
		cell.cellColor[1] += random_color_val;

		random_val = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
		random_color_val = 0 + random_val * (1 - 0);
		cell.cellColor[2] += random_color_val;
        voronoiCells.push_back(cell);
    }

    return voronoiCells;
}


// ----------------------------------------------------------------------------------
// RENDER CODE. If you need to edit any of the code here, feel free to do so.
// ----------------------------------------------------------------------------------

#include <cstddef>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <stdexcept>


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
const int WINDOW_WIDTH = 1280;
const int WINDOW_HEIGHT = 720;
struct AppData {
	std::vector<Sites> testCases;
	int currentTestCase;
	GLuint vbo;
	size_t numVertices;
	bool showTriangulation;
	Point player;
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

	//EDIT THIS CODE IF YOU WANT EXTRA KEYBOARD INPUTS

	AppData* appData = reinterpret_cast<AppData*>(glfwGetWindowUserPointer(window));
	if (action == GLFW_PRESS) {
		if (key == GLFW_KEY_LEFT) {
			appData->currentTestCase = ((appData->currentTestCase - 1) + appData->testCases.size()) % appData->testCases.size();
		}
		else if (key == GLFW_KEY_RIGHT) {
			appData->currentTestCase = (appData->currentTestCase + 1) % appData->testCases.size();
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

// copy pasted from previous exercise
void AppendShapeVertices(const std::vector<Point>& points, std::vector<Vertex>& vertices, float zOrder = 0.0f, float r = 1.0f, float g = 1.0f, float b = 1.0f) {
	if (points.size() < 2) {
		return;
	}
	for (size_t i = 2; i < points.size(); ++i) {
		vertices.push_back( { points[0].x, points[0].y, zOrder, r, g, b } );
		vertices.push_back( { points[i - 1].x, points[i - 1].y, zOrder, r, g, b } );
		vertices.push_back( { points[i].x, points[i].y, zOrder, r, g, b } );
	}
}


// The RefreshScene() function draws the sites and cells of the Voronoi, so if you want to
// make edits to the render, edit that function.

void RefreshScene(AppData* data) {

	//EDIT THIS CODE FOR ANY NEW RENDER

	std::vector<Vertex> vertices;
	// Draw sites as a circle
	Sites& sites = data->testCases[data->currentTestCase];
	for (size_t i = 0; i < sites.sitelist.size(); ++i) {
		AppendCircleVertices(sites.sitelist[i], 5.0f, vertices, 0.0f, 1.0f, 1.0f, 1.0f);
	}

	// Get the data of the voronoi cells
	std::vector<Cell> voronoiCells = VoronoiDiagram(data->testCases[data->currentTestCase].sitelist);

	
	for (Cell cell : voronoiCells) {

        for (int i = 0; i < cell.vertices.size(); ++i) {
            Point p0 = cell.vertices[i];
            Point p1 = cell.vertices[(i + 1) % cell.vertices.size()];
            AppendLineVertices(p0, p1, 2.0f, vertices, 0.0f, cell.edgeColor[0], cell.edgeColor[1], cell.edgeColor[2]);
			// add vertex polygon to shade
			AppendShapeVertices(cell.vertices, vertices, 0.0f, cell.cellColor[0], cell.cellColor[1], cell.cellColor[2]);
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
	appData.player = { WINDOW_WIDTH / 2.0f, WINDOW_HEIGHT / 2.0f };
	std::ifstream file(INPUT_FILE);
	if (file.fail())
	{
		std::cout << "Failed to read file test.txt" << std::endl;
		return 1;
	}

	size_t numTestCases;
	file >> numTestCases;
	for (size_t testCaseNumber = 0; testCaseNumber < numTestCases; ++testCaseNumber) {
		Sites sites;
		
		size_t numPoints;
		file >> numPoints;

		std::vector<Point> points;
		for (size_t j = 0; j < numPoints; ++j) {
			Point point;
			file >> point.x >> point.y;
			points.push_back(point);
		}
		sites.sitelist = points;

		appData.testCases.push_back(sites);
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

	GLFWwindow* window = glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "LastNameLastNameLastName - Voronoi Diagram", nullptr, nullptr);
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
