# Double Pendulum
A system of two rods and masses, with angles (θ1,θ2). Lagrangian mechanics yields coupled ODEs, solved numerically. Tiny changes in initial angles produce huge variations. Adjust code for mesmerizing chaos!

## The Math Explained

We label:
- θ₁ as the upper bob’s angle from vertical
- θ₂ as the lower bob’s angle from vertical
- ω₁ = dθ₁/dt and ω₂ = dθ₂/dt
- L₁, L₂ as rod lengths
- m₁, m₂ as bob masses
- g as gravitational acceleration

### 1. Kinetic Energy (T)
Positions:
- x₁ = L₁ sin(θ₁), y₁ = −L₁ cos(θ₁)
- x₂ = x₁ + L₂ sin(θ₂), y₂ = y₁ − L₂ cos(θ₂)

Velocities:
- ẋ₁ = L₁ cos(θ₁) ω₁,  ẏ₁ = L₁ sin(θ₁) ω₁
- ẋ₂ = ẋ₁ + L₂ cos(θ₂) ω₂,  ẏ₂ = ẏ₁ + L₂ sin(θ₂) ω₂

Hence:
T = ½m₁(ẋ₁² + ẏ₁²) + ½m₂(ẋ₂² + ẏ₂²)

### 2. Potential Energy (V)
Choose the pivot as reference:
V = m₁g(L₁ − y₁) + m₂g(L₁ + L₂ − y₂)

### 3. Equations of Motion
From L = T − V and Euler-Lagrange,
dθ₁/dt = ω₁,  
dθ₂/dt = ω₂,

dω₁/dt = 
[−g(2m₁ + m₂) sin(θ₁) − m₂g sin(θ₁ − 2θ₂) 
− 2m₂ sin(θ₁ − θ₂)(L₂ω₂² + L₁ω₁² cos(θ₁ − θ₂))]
/ [L₁(2m₁ + m₂ − m₂ cos(2θ₁ − 2θ₂))],

dω₂/dt = 
[2 sin(θ₁ − θ₂)(L₁ω₁²(m₁ + m₂) + g(m₁ + m₂) cos(θ₁) + L₂m₂ω₂² cos(θ₁ − θ₂))]
/ [L₂(2m₁ + m₂ − m₂ cos(2θ₁ − 2θ₂))].

Solving these ODEs numerically reveals the pendulum’s chaotic motion.
