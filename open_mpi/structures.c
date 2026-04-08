// Uniform charge density

typedef struct {
    float epsilon; // divide by 10 power -12
    float uniform_charge_density; // in pC per m
    
} Physics;

typedef struct {
    float step_size;
    float iteration;
} Simulation;

typedef struct {
    float left;
    float right;

    float up;
    float down;

    float top;
    float bottom;
} Boundary;

typedef struct {
    int x;
    int y;
    int z;
} Cartesian_int;

typedef struct {
    float x;
    float y;
    float z;
} Cartesian_float;

typedef struct {
        float* vertical_1;
        float* vertical_2;
        float* horizontal_1;
        float* horizontal_2;
    } Buffer;