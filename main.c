#include <raylib.h>
#include <raymath.h>
#include <math.h>
#include <stdio.h>

#define G 200.0f // Gravitational constant
#define C_SPEED 1000.0f // Speed of light approximation
#define MAX_BODIES 500
#define TRAIL_LENGTH 200
#define DISK_INNER_RADIUS 80.0f
#define DISK_OUTER_RADIUS 600.0f

typedef struct {
    Vector2 position;
    Vector2 velocity;
    float mass;
    float radius;
    Color color;
    Vector2 trail[TRAIL_LENGTH];
    int trailIndex;
    bool active;
} Body;

typedef struct {
    Vector2 position;
    Vector2 velocity;
} State;

typedef struct {
    Vector2 dPosition; // velocity
    Vector2 dVelocity; // acceleration
} Derivative;

// Extract state from bodies for RK4
void getStates(Body bodies[], State states[]) {
    for(int i=0; i<MAX_BODIES; i++) {
        states[i].position = bodies[i].position;
        states[i].velocity = bodies[i].velocity;
    }
}

// Calculate derivatives (accelerations) for a given state
void calculateDerivatives(State states[], Derivative derivs[], Body bodies[], bool enableDrag) {
    for (int i = 0; i < MAX_BODIES; i++) {
        if (!bodies[i].active) {
            derivs[i].dPosition = (Vector2){0,0};
            derivs[i].dVelocity = (Vector2){0,0};
            continue;
        }

        derivs[i].dPosition = states[i].velocity; // dr/dt = v
        Vector2 force = { 0.0f, 0.0f };

        // Drag Force (Gas/Dust) - Accretion Disk Model
        if (enableDrag && bodies[0].active) {
            // Assume Body 0 is the center (Sun)
            float distToCenter = Vector2Distance(states[i].position, states[0].position);
            
            if (distToCenter > DISK_INNER_RADIUS && distToCenter < DISK_OUTER_RADIUS) {
                // Density function: Higher density closer to center
                // rho(r) = 1 - (r - Rin) / (Rout - Rin)
                float normalizedDist = (distToCenter - DISK_INNER_RADIUS) / (DISK_OUTER_RADIUS - DISK_INNER_RADIUS);
                float density = 1.0f - normalizedDist;
                if (density < 0.0f) density = 0.0f;
                
                // Reduced drag coefficient significantly
                // Was 0.05f, now 0.002f to be much more subtle
                float dragCoeff = 0.002f * density; 
                
                // F_drag = -c * rho * v
                Vector2 drag = Vector2Scale(states[i].velocity, -dragCoeff);
                force = Vector2Add(force, drag);
            }
        }

        for (int j = 0; j < MAX_BODIES; j++) {
            if (i == j || !bodies[j].active) continue;

            Vector2 direction = Vector2Subtract(states[j].position, states[i].position);
            float distance = Vector2Length(direction);
            
            // Softening to avoid singularities, though collision should handle this
            float minDist = bodies[i].radius + bodies[j].radius; 
            if (distance < minDist) distance = minDist; 

            float forceMagnitude = (G * bodies[i].mass * bodies[j].mass) / (distance * distance);

            // Relativistic correction (Precession)
            // F_gr = F_newton * (1 + 3(h/rc)^2)
            Vector2 relVel = Vector2Subtract(states[j].velocity, states[i].velocity);
            float h = direction.x * relVel.y - direction.y * relVel.x; // 2D cross product magnitude
            float correction = (3.0f * h * h) / (C_SPEED * C_SPEED * distance * distance);
            forceMagnitude *= (1.0f + correction);

            Vector2 forceVec = Vector2Scale(Vector2Normalize(direction), forceMagnitude);
            force = Vector2Add(force, forceVec);
        }
        derivs[i].dVelocity = Vector2Scale(force, 1.0f / bodies[i].mass); // dv/dt = a
    }
}

// Runge-Kutta 4th Order Integration
void integrateRK4(Body bodies[], float dt, bool enableDrag) {
    State initialStates[MAX_BODIES];
    getStates(bodies, initialStates);

    Derivative k1[MAX_BODIES], k2[MAX_BODIES], k3[MAX_BODIES], k4[MAX_BODIES];
    State tempStates[MAX_BODIES];

    // k1
    calculateDerivatives(initialStates, k1, bodies, enableDrag);

    // k2
    for(int i=0; i<MAX_BODIES; i++) {
        tempStates[i].position = Vector2Add(initialStates[i].position, Vector2Scale(k1[i].dPosition, dt * 0.5f));
        tempStates[i].velocity = Vector2Add(initialStates[i].velocity, Vector2Scale(k1[i].dVelocity, dt * 0.5f));
    }
    calculateDerivatives(tempStates, k2, bodies, enableDrag);

    // k3
    for(int i=0; i<MAX_BODIES; i++) {
        tempStates[i].position = Vector2Add(initialStates[i].position, Vector2Scale(k2[i].dPosition, dt * 0.5f));
        tempStates[i].velocity = Vector2Add(initialStates[i].velocity, Vector2Scale(k2[i].dVelocity, dt * 0.5f));
    }
    calculateDerivatives(tempStates, k3, bodies, enableDrag);

    // k4
    for(int i=0; i<MAX_BODIES; i++) {
        tempStates[i].position = Vector2Add(initialStates[i].position, Vector2Scale(k3[i].dPosition, dt));
        tempStates[i].velocity = Vector2Add(initialStates[i].velocity, Vector2Scale(k3[i].dVelocity, dt));
    }
    calculateDerivatives(tempStates, k4, bodies, enableDrag);

    // Update
    for(int i=0; i<MAX_BODIES; i++) {
        if (!bodies[i].active) continue;
        
        Vector2 dPos = Vector2Scale(Vector2Add(Vector2Add(k1[i].dPosition, Vector2Scale(k2[i].dPosition, 2.0f)), Vector2Add(Vector2Scale(k3[i].dPosition, 2.0f), k4[i].dPosition)), dt / 6.0f);
        Vector2 dVel = Vector2Scale(Vector2Add(Vector2Add(k1[i].dVelocity, Vector2Scale(k2[i].dVelocity, 2.0f)), Vector2Add(Vector2Scale(k3[i].dVelocity, 2.0f), k4[i].dVelocity)), dt / 6.0f);

        bodies[i].position = Vector2Add(bodies[i].position, dPos);
        bodies[i].velocity = Vector2Add(bodies[i].velocity, dVel);
    }
}

void handleCollisions(Body bodies[]) {
    for (int i = 0; i < MAX_BODIES; i++) {
        if (!bodies[i].active) continue;
        for (int j = i + 1; j < MAX_BODIES; j++) {
            if (!bodies[j].active) continue;

            if (CheckCollisionCircles(bodies[i].position, bodies[i].radius, bodies[j].position, bodies[j].radius)) {
                // Merge j into i (Fusion)
                Body *b1 = &bodies[i];
                Body *b2 = &bodies[j];

                // Conservation of momentum: (m1v1 + m2v2) / (m1 + m2)
                Vector2 momentum1 = Vector2Scale(b1->velocity, b1->mass);
                Vector2 momentum2 = Vector2Scale(b2->velocity, b2->mass);
                Vector2 totalMomentum = Vector2Add(momentum1, momentum2);
                float totalMass = b1->mass + b2->mass;

                b1->velocity = Vector2Scale(totalMomentum, 1.0f / totalMass);
                
                // Weighted position
                b1->position = Vector2Scale(Vector2Add(Vector2Scale(b1->position, b1->mass), Vector2Scale(b2->position, b2->mass)), 1.0f/totalMass);

                // New radius (Volume conservation: r^3 = r1^3 + r2^3)
                // This is much more realistic for 3D bodies and prevents rapid "2D area" expansion
                // which can cause chain reactions.
                b1->radius = cbrtf(powf(b1->radius, 3.0f) + powf(b2->radius, 3.0f));
                b1->mass = totalMass;
                
                // Deactivate b2
                b2->active = false;
            }
        }
    }
}

void explodeBody(Body bodies[], int index) {
    // Deactivate original
    bodies[index].active = false;
    
    int fragments = 8;
    float newMass = bodies[index].mass / fragments;
    float newRadius = bodies[index].radius / 2.0f;
    if (newRadius < 2.0f) newRadius = 2.0f;
    
    for (int k = 0; k < fragments; k++) {
        // Find free slot
        for (int j = 0; j < MAX_BODIES; j++) {
            if (!bodies[j].active) {
                bodies[j].active = true;
                bodies[j].mass = newMass;
                bodies[j].radius = newRadius;
                // Vary color slightly
                bodies[j].color = (Color){ 
                    (unsigned char)fminf(255, bodies[index].color.r + GetRandomValue(-20, 20)),
                    (unsigned char)fminf(255, bodies[index].color.g + GetRandomValue(-20, 20)),
                    (unsigned char)fminf(255, bodies[index].color.b + GetRandomValue(-20, 20)),
                    255 
                };
                
                // Random offset
                Vector2 offset = { (float)GetRandomValue(-5, 5), (float)GetRandomValue(-5, 5) };
                bodies[j].position = Vector2Add(bodies[index].position, offset);
                
                // Random velocity spread
                Vector2 velSpread = { (float)GetRandomValue(-20, 20) / 10.0f, (float)GetRandomValue(-20, 20) / 10.0f };
                bodies[j].velocity = Vector2Add(bodies[index].velocity, velSpread);
                
                // Reset trail
                bodies[j].trailIndex = 0;
                for(int t=0; t<TRAIL_LENGTH; t++) bodies[j].trail[t] = bodies[j].position;
                break;
            }
        }
    }
}

void checkRocheLimit(Body bodies[]) {
    // Check all bodies against Sun (body 0)
    if (!bodies[0].active) return;
    
    for (int i = 1; i < MAX_BODIES; i++) {
        if (!bodies[i].active) continue;
        
        float dist = Vector2Distance(bodies[i].position, bodies[0].position);
        
        // Roche limit approximation: d = 2.44 * R_sat * (M_sun / M_sat)^(1/3)
        // Avoid division by zero
        if (bodies[i].mass < 0.1f) continue; 
        
        float massRatio = bodies[0].mass / bodies[i].mass;
        // Use Rigid Body Roche Limit (1.26) instead of Fluid (2.44)
        // This prevents planets from exploding too early/far from Sun
        float rocheLimit = 1.26f * bodies[i].radius * cbrtf(massRatio);
        
        // Add a safety margin or ensure we don't explode things that are already tiny fragments
        // Let's say fragments with radius < 3.0 don't explode further to prevent infinite recursion/lag
        if (bodies[i].radius > 3.0f && dist < rocheLimit) {
            explodeBody(bodies, i);
        }
    }
}

void drawLagrangePoints(Body bodies[]) {
    // Assuming Body 0 is Sun. Find the most massive planet.
    if (!bodies[0].active) return;

    int heaviestIndex = -1;
    float maxMass = 0.0f;

    for (int i = 1; i < MAX_BODIES; i++) {
        if (bodies[i].active && bodies[i].mass > maxMass) {
            maxMass = bodies[i].mass;
            heaviestIndex = i;
        }
    }

    if (heaviestIndex == -1) return;

    Body *m1 = &bodies[0];
    Body *m2 = &bodies[heaviestIndex];

    Vector2 r1 = m1->position;
    Vector2 r2 = m2->position;
    Vector2 R_vec = Vector2Subtract(r2, r1);
    float R = Vector2Length(R_vec);
    Vector2 u = Vector2Scale(R_vec, 1.0f / R);

    float massRatio = m2->mass / m1->mass; // q = m2/m1
    // Hill Radius approximation: r_H = R * (m2 / 3m1)^(1/3)
    float hillRadius = R * cbrtf(massRatio / 3.0f);

    // L1: Between M1 and M2, close to M2
    Vector2 l1 = Vector2Subtract(r2, Vector2Scale(u, hillRadius));
    
    // L2: Beyond M2
    Vector2 l2 = Vector2Add(r2, Vector2Scale(u, hillRadius));

    // L3: Opposite side of M1
    // Approx: R * (1 + 5/12 * m2/m1)
    float l3_dist = R * (1.0f + (5.0f/12.0f) * massRatio);
    Vector2 l3 = Vector2Subtract(r1, Vector2Scale(u, l3_dist));

    // L4 & L5: Equilateral triangle with M1 and M2
    // Rotate R_vec by 60 degrees (+/- PI/3)
    
    // Manual rotation since Vector2Rotate might vary in availability/version
    float cos60 = 0.5f;
    float sin60 = 0.8660254f;

    Vector2 l4_offset = {
        R_vec.x * cos60 - R_vec.y * sin60,
        R_vec.x * sin60 + R_vec.y * cos60
    };
    Vector2 l4 = Vector2Add(r1, l4_offset);

    Vector2 l5_offset = {
        R_vec.x * cos60 - R_vec.y * (-sin60),
        R_vec.x * (-sin60) + R_vec.y * cos60
    };
    Vector2 l5 = Vector2Add(r1, l5_offset);

    // Draw
    Color lColor = VIOLET;
    float lRadius = 5.0f;
    
    DrawCircleV(l1, lRadius, lColor); DrawText("L1", l1.x + 5, l1.y - 10, 10, WHITE);
    DrawCircleV(l2, lRadius, lColor); DrawText("L2", l2.x + 5, l2.y - 10, 10, WHITE);
    DrawCircleV(l3, lRadius, lColor); DrawText("L3", l3.x + 5, l3.y - 10, 10, WHITE);
    DrawCircleV(l4, lRadius, lColor); DrawText("L4", l4.x + 5, l4.y - 10, 10, WHITE);
    DrawCircleV(l5, lRadius, lColor); DrawText("L5", l5.x + 5, l5.y - 10, 10, WHITE);
}

void initBodies(Body bodies[], int screenWidth, int screenHeight) {
    // Clear all bodies
    for(int i=0; i<MAX_BODIES; i++) bodies[i].active = false;

    // Sun
    bodies[0] = (Body){
        .position = { screenWidth / 2.0f, screenHeight / 2.0f },
        .velocity = { 0.0f, 0.0f },
        .mass = 10000.0f,
        .radius = 40.0f,
        .color = YELLOW,
        .active = true
    };
    // Planet 1 (Earth-like)
    bodies[1] = (Body){
        .position = { screenWidth / 2.0f + 200.0f, screenHeight / 2.0f },
        .velocity = { 0.0f, 100.0f },
        .mass = 10.0f,
        .radius = 10.0f,
        .color = BLUE,
        .active = true
    };
    // Planet 2 (Mars-like)
    bodies[2] = (Body){
        .position = { screenWidth / 2.0f + 300.0f, screenHeight / 2.0f },
        .velocity = { 0.0f, 81.6f },
        .mass = 8.0f,
        .radius = 8.0f,
        .color = RED,
        .active = true
    };
    // Planet 3 (Jupiter-like)
    bodies[3] = (Body){
        .position = { screenWidth / 2.0f + 450.0f, screenHeight / 2.0f },
        .velocity = { 0.0f, 66.6f },
        .mass = 50.0f,
        .radius = 20.0f,
        .color = ORANGE,
        .active = true
    };

    // Initialize trails
    for (int i = 0; i < MAX_BODIES; i++) {
        for (int j = 0; j < TRAIL_LENGTH; j++) {
            bodies[i].trail[j] = bodies[i].position;
        }
        bodies[i].trailIndex = 0;
    }
}

typedef enum {
    STATE_SIMULATION,
    STATE_MENU,
    STATE_SETTINGS
} AppState;

int main(void)
{
    const int screenWidth = 1200;
    const int screenHeight = 800;

    InitWindow(screenWidth, screenHeight, "Solar System Simulation - Solaray");
    SetTargetFPS(60);
    SetExitKey(KEY_NULL); // Disable default ESC behavior

    // Camera setup
    Camera2D camera = { 0 };
    camera.zoom = 1.0f;
    camera.target = (Vector2){ screenWidth / 2.0f, screenHeight / 2.0f };
    camera.offset = (Vector2){ screenWidth / 2.0f, screenHeight / 2.0f };
    camera.rotation = 0.0f;

    float timeScale = 1.0f;
    bool enableDrag = false;
    bool showLagrange = false;
    bool enableRoche = true; // New setting

    Body bodies[MAX_BODIES];
    initBodies(bodies, screenWidth, screenHeight);

    // Creation Mode State
    bool creationMode = false;
    bool isDragging = false;
    Vector2 dragStartPos = {0};
    float newBodyMass = 10.0f;

    AppState currentState = STATE_SIMULATION;
    bool shouldExit = false;

    while (!shouldExit && !WindowShouldClose()) {
        // Global Input
        if (IsKeyPressed(KEY_ESCAPE)) {
            if (currentState == STATE_SIMULATION) currentState = STATE_MENU;
            else if (currentState == STATE_MENU) currentState = STATE_SIMULATION;
            else if (currentState == STATE_SETTINGS) currentState = STATE_MENU;
        }

        if (currentState == STATE_SIMULATION) {
            // Input Handling
            if (IsKeyPressed(KEY_N)) creationMode = !creationMode;

            if (creationMode) {
                // Creation Logic
                Vector2 mouseWorldPos = GetScreenToWorld2D(GetMousePosition(), camera);
                
                // Adjust Mass with Scroll
                float wheel = GetMouseWheelMove();
                if (wheel != 0) {
                    newBodyMass += wheel * 2.0f;
                    if (newBodyMass < 1.0f) newBodyMass = 1.0f;
                }

                if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) {
                    isDragging = true;
                    dragStartPos = mouseWorldPos;
                }

                if (IsMouseButtonReleased(MOUSE_BUTTON_LEFT) && isDragging) {
                    isDragging = false;
                    // Spawn Body
                    for (int i = 0; i < MAX_BODIES; i++) {
                        if (!bodies[i].active) {
                            bodies[i].active = true;
                            bodies[i].position = dragStartPos;
                            bodies[i].velocity = Vector2Subtract(mouseWorldPos, dragStartPos);
                            bodies[i].mass = newBodyMass;
                            bodies[i].radius = sqrtf(newBodyMass) * 3.0f; // Rough visual scaling
                            if (bodies[i].radius < 5.0f) bodies[i].radius = 5.0f;
                            bodies[i].color = (Color){ GetRandomValue(100, 255), GetRandomValue(100, 255), GetRandomValue(100, 255), 255 };
                            
                            // Reset trail
                            for(int t=0; t<TRAIL_LENGTH; t++) bodies[i].trail[t] = bodies[i].position;
                            bodies[i].trailIndex = 0;
                            break;
                        }
                    }
                }
            } else {
                // Normal Camera Controls
                float wheel = GetMouseWheelMove();
                if (wheel != 0) {
                    Vector2 mouseWorldPos = GetScreenToWorld2D(GetMousePosition(), camera);
                    camera.offset = GetMousePosition();
                    camera.target = mouseWorldPos;
                    float scaleFactor = 1.0f + (0.25f * fabs(wheel));
                    if (wheel < 0) camera.zoom = 1.0f / scaleFactor * camera.zoom;
                    else camera.zoom = scaleFactor * camera.zoom;
                    if (camera.zoom < 0.1f) camera.zoom = 0.1f;
                }
                if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT) || IsMouseButtonDown(MOUSE_BUTTON_MIDDLE)) {
                    Vector2 delta = GetMouseDelta();
                    delta = Vector2Scale(delta, -1.0f / camera.zoom);
                    camera.target = Vector2Add(camera.target, delta);
                }
            }

            if (IsKeyPressed(KEY_RIGHT)) timeScale *= 2.0f;
            if (IsKeyPressed(KEY_LEFT)) timeScale *= 0.5f;
            if (IsKeyPressed(KEY_SPACE)) timeScale = (timeScale == 0.0f) ? 1.0f : 0.0f;

            float dt = GetFrameTime() * timeScale;
            
            // Physics Loop
            int subSteps = (int)ceilf(fabs(timeScale)); 
            if (subSteps < 1) subSteps = 1;
            float subDt = dt / subSteps;

            for (int step = 0; step < subSteps; step++) {
                integrateRK4(bodies, subDt, enableDrag);
                handleCollisions(bodies);
                if (enableRoche) checkRocheLimit(bodies);
            }

            // Update Trails
            if (timeScale != 0.0f) {
                for (int i = 0; i < MAX_BODIES; i++) {
                    if (!bodies[i].active) continue;
                    bodies[i].trail[bodies[i].trailIndex] = bodies[i].position;
                    bodies[i].trailIndex = (bodies[i].trailIndex + 1) % TRAIL_LENGTH;
                }
            }
        }

        // Rendering
        BeginDrawing();
        ClearBackground(BLACK);
        
        BeginMode2D(camera);
        
        // Draw Accretion Disk
        if (enableDrag && bodies[0].active) {
            // Draw a series of rings to simulate the disk gradient
            int rings = 20;
            float step = (DISK_OUTER_RADIUS - DISK_INNER_RADIUS) / rings;
            for (int i = 0; i < rings; i++) {
                float rInner = DISK_INNER_RADIUS + i * step;
                float rOuter = rInner + step;
                float alpha = 0.3f * (1.0f - (float)i / rings); // Fade out
                Color ringColor = (Color){ 100, 100, 255, (unsigned char)(alpha * 255) };
                DrawRing(bodies[0].position, rInner, rOuter, 0, 360, 64, ringColor);
            }
        }

        // Draw Trails
        for (int i = 0; i < MAX_BODIES; i++) {
            if (!bodies[i].active) continue;
            for (int j = 0; j < TRAIL_LENGTH - 1; j++) {
                int idx = (bodies[i].trailIndex + j) % TRAIL_LENGTH;
                int nextIdx = (idx + 1) % TRAIL_LENGTH;
                if (Vector2DistanceSqr(bodies[i].trail[idx], bodies[i].trail[nextIdx]) > 0.001f)
                    DrawLineV(bodies[i].trail[idx], bodies[i].trail[nextIdx], Fade(bodies[i].color, 0.4f));
            }
        }
        
        // Draw Bodies
        for (int i = 0; i < MAX_BODIES; i++) {
            if (!bodies[i].active) continue;
            DrawCircleV(bodies[i].position, bodies[i].radius, bodies[i].color);
        }

        // Creation Mode Visualization
        if (creationMode && currentState == STATE_SIMULATION) {
            Vector2 mouseWorldPos = GetScreenToWorld2D(GetMousePosition(), camera);
            if (isDragging) {
                DrawLineV(dragStartPos, mouseWorldPos, WHITE);
                DrawCircleV(dragStartPos, 5.0f, GREEN); // Start point
                // Preview Body
                float previewRadius = sqrtf(newBodyMass) * 3.0f;
                if (previewRadius < 5.0f) previewRadius = 5.0f;
                DrawCircleLinesV(dragStartPos, previewRadius, Fade(WHITE, 0.5f));
                
                // Velocity Text
                Vector2 vel = Vector2Subtract(mouseWorldPos, dragStartPos);
                DrawText(TextFormat("Vel: %.1f", Vector2Length(vel)), mouseWorldPos.x + 10, mouseWorldPos.y, 10, WHITE);
            } else {
                // Cursor Preview
                float previewRadius = sqrtf(newBodyMass) * 3.0f;
                if (previewRadius < 5.0f) previewRadius = 5.0f;
                DrawCircleLinesV(mouseWorldPos, previewRadius, Fade(GRAY, 0.5f));
            }
        }

        if (showLagrange) {
            drawLagrangePoints(bodies);
        }
        
        EndMode2D();
        
        // Hover Info (Screen Space)
        if (currentState == STATE_SIMULATION) {
            Vector2 mousePos = GetMousePosition();
            Vector2 mouseWorldPos = GetScreenToWorld2D(mousePos, camera);
            for (int i = 0; i < MAX_BODIES; i++) {
                if (!bodies[i].active) continue;
                if (CheckCollisionPointCircle(mouseWorldPos, bodies[i].position, bodies[i].radius)) {
                    float speed = Vector2Length(bodies[i].velocity);
                    float distToSun = Vector2Distance(bodies[i].position, bodies[0].position);
                    float kineticE = 0.5f * bodies[i].mass * speed * speed;
                    
                    char infoText[512];
                    sprintf(infoText, "Mass: %.1f\nSpeed: %.1f\nDist to Sun: %.1f\nKinetic E: %.1e\nPos: (%.0f, %.0f)", 
                            bodies[i].mass, speed, distToSun, kineticE, bodies[i].position.x, bodies[i].position.y);
                    
                    Vector2 screenPos = GetWorldToScreen2D(bodies[i].position, camera);
                    float screenRadius = bodies[i].radius * camera.zoom;

                    DrawRectangle(screenPos.x + screenRadius + 5, screenPos.y - 60, 200, 100, Fade(DARKGRAY, 0.9f));
                    DrawRectangleLines(screenPos.x + screenRadius + 5, screenPos.y - 60, 200, 100, WHITE);
                    DrawText(infoText, screenPos.x + screenRadius + 10, screenPos.y - 55, 10, WHITE);
                }
            }
        }
        
        // Global Stats Calculation
        float totalEnergy = 0.0f;
        float totalMomentum = 0.0f;
        int activeCount = 0;
        for(int i=0; i<MAX_BODIES; i++) {
            if(!bodies[i].active) continue;
            activeCount++;
            float v = Vector2Length(bodies[i].velocity);
            totalEnergy += 0.5f * bodies[i].mass * v * v; // Kinetic
            totalMomentum += bodies[i].mass * v;
            
            // Potential (simplified pair-wise)
            for(int j=i+1; j<MAX_BODIES; j++) {
                if(!bodies[j].active) continue;
                float r = Vector2Distance(bodies[i].position, bodies[j].position);
                if (r > 1.0f) totalEnergy -= (G * bodies[i].mass * bodies[j].mass) / r;
            }
        }

        DrawFPS(10, 10);
        DrawText(TextFormat("Time Scale: %.2fx", timeScale), 10, 30, 20, WHITE);
        
        // Static Features
        DrawText("RK4 | N-Body | Collisions | Relativistic Precession", 10, 50, 20, GREEN);
        
        // Dynamic Features Status
        int statusY = 70;
        int x = 10;
        
        DrawText("Roche Limit:", x, statusY, 20, LIGHTGRAY);
        x += MeasureText("Roche Limit: ", 20);
        DrawText(enableRoche ? "ON" : "OFF", x, statusY, 20, enableRoche ? GREEN : RED);
        x += MeasureText("ON ", 20) + 10;
        
        DrawText("| Accretion:", x, statusY, 20, LIGHTGRAY);
        x += MeasureText("| Accretion: ", 20);
        DrawText(enableDrag ? "ON" : "OFF", x, statusY, 20, enableDrag ? GREEN : RED);
        x += MeasureText("ON ", 20) + 10;

        DrawText("| Lagrange:", x, statusY, 20, LIGHTGRAY);
        x += MeasureText("| Lagrange: ", 20);
        DrawText(showLagrange ? "ON" : "OFF", x, statusY, 20, showLagrange ? GREEN : RED);
        x += MeasureText("ON ", 20) + 10;

        DrawText("| Create (N):", x, statusY, 20, LIGHTGRAY);
        x += MeasureText("| Create (N): ", 20);
        DrawText(creationMode ? "ON" : "OFF", x, statusY, 20, creationMode ? GREEN : RED);
        
        if (creationMode && currentState == STATE_SIMULATION) {
            DrawText(TextFormat("Creation Mode: Click & Drag to launch. Scroll to set Mass: %.1f", newBodyMass), 10, 100, 20, YELLOW);
        }

        // Global Stats Display
        int statsY = screenHeight - 120;
        DrawRectangle(0, statsY, 300, 120, Fade(BLACK, 0.5f));
        DrawText("Global Statistics", 10, statsY + 10, 20, SKYBLUE);
        DrawText(TextFormat("Active Bodies: %d", activeCount), 10, statsY + 40, 10, WHITE);
        DrawText(TextFormat("Total Energy: %.2e", totalEnergy), 10, statsY + 60, 10, WHITE);
        DrawText(TextFormat("Total Momentum: %.2e", totalMomentum), 10, statsY + 80, 10, WHITE);

        // Menu Overlay
        if (currentState == STATE_MENU) {
            DrawRectangle(0, 0, screenWidth, screenHeight, Fade(BLACK, 0.7f));
            int menuX = screenWidth / 2 - 100;
            int menuY = screenHeight / 2 - 100;
            
            DrawText("PAUSED", screenWidth/2 - MeasureText("PAUSED", 40)/2, menuY - 60, 40, WHITE);

            // Simple Button Logic
            Vector2 mouse = GetMousePosition();
            
            // Resume
            Rectangle btnResume = { menuX, menuY, 200, 40 };
            bool hoverResume = CheckCollisionPointRec(mouse, btnResume);
            DrawRectangleRec(btnResume, hoverResume ? GRAY : DARKGRAY);
            DrawText("Resume", btnResume.x + 20, btnResume.y + 10, 20, WHITE);
            if (hoverResume && IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) currentState = STATE_SIMULATION;

            // Reset
            Rectangle btnReset = { menuX, menuY + 50, 200, 40 };
            bool hoverReset = CheckCollisionPointRec(mouse, btnReset);
            DrawRectangleRec(btnReset, hoverReset ? GRAY : DARKGRAY);
            DrawText("Reset", btnReset.x + 20, btnReset.y + 10, 20, WHITE);
            if (hoverReset && IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) {
                initBodies(bodies, screenWidth, screenHeight);
                currentState = STATE_SIMULATION;
            }

            // Settings
            Rectangle btnSettings = { menuX, menuY + 100, 200, 40 };
            bool hoverSettings = CheckCollisionPointRec(mouse, btnSettings);
            DrawRectangleRec(btnSettings, hoverSettings ? GRAY : DARKGRAY);
            DrawText("Settings", btnSettings.x + 20, btnSettings.y + 10, 20, WHITE);
            if (hoverSettings && IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) currentState = STATE_SETTINGS;

            // Quit
            Rectangle btnQuit = { menuX, menuY + 150, 200, 40 };
            bool hoverQuit = CheckCollisionPointRec(mouse, btnQuit);
            DrawRectangleRec(btnQuit, hoverQuit ? RED : MAROON);
            DrawText("Quit", btnQuit.x + 20, btnQuit.y + 10, 20, WHITE);
            if (hoverQuit && IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) shouldExit = true;
        }

        // Settings Overlay
        else if (currentState == STATE_SETTINGS) {
            DrawRectangle(0, 0, screenWidth, screenHeight, Fade(BLACK, 0.8f));
            int menuX = screenWidth / 2 - 150;
            int menuY = screenHeight / 2 - 100;
            
            DrawText("SETTINGS", screenWidth/2 - MeasureText("SETTINGS", 40)/2, menuY - 60, 40, WHITE);
            Vector2 mouse = GetMousePosition();

            // Toggle Roche
            Rectangle btnRoche = { menuX, menuY, 300, 40 };
            bool hoverRoche = CheckCollisionPointRec(mouse, btnRoche);
            DrawRectangleRec(btnRoche, hoverRoche ? GRAY : DARKGRAY);
            DrawText(TextFormat("Roche Limit: %s", enableRoche ? "ON" : "OFF"), btnRoche.x + 20, btnRoche.y + 10, 20, enableRoche ? GREEN : RED);
            if (hoverRoche && IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) enableRoche = !enableRoche;

            // Toggle Drag
            Rectangle btnDrag = { menuX, menuY + 50, 300, 40 };
            bool hoverDrag = CheckCollisionPointRec(mouse, btnDrag);
            DrawRectangleRec(btnDrag, hoverDrag ? GRAY : DARKGRAY);
            DrawText(TextFormat("Accretion Drag: %s", enableDrag ? "ON" : "OFF"), btnDrag.x + 20, btnDrag.y + 10, 20, enableDrag ? GREEN : RED);
            if (hoverDrag && IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) enableDrag = !enableDrag;

            // Toggle Lagrange
            Rectangle btnLag = { menuX, menuY + 100, 300, 40 };
            bool hoverLag = CheckCollisionPointRec(mouse, btnLag);
            DrawRectangleRec(btnLag, hoverLag ? GRAY : DARKGRAY);
            DrawText(TextFormat("Lagrange Points: %s", showLagrange ? "ON" : "OFF"), btnLag.x + 20, btnLag.y + 10, 20, showLagrange ? GREEN : RED);
            if (hoverLag && IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) showLagrange = !showLagrange;

            // Back
            Rectangle btnBack = { menuX, menuY + 160, 300, 40 };
            bool hoverBack = CheckCollisionPointRec(mouse, btnBack);
            DrawRectangleRec(btnBack, hoverBack ? GRAY : DARKGRAY);
            DrawText("Back", btnBack.x + 120, btnBack.y + 10, 20, WHITE);
            if (hoverBack && IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) currentState = STATE_MENU;
        }

        EndDrawing();
    }
    CloseWindow();
    return 0;
}
