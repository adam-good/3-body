# Three Body Problem as Defined By https://en.wikipedia.org/wiki/Three-body_problem
module ThreeBody

export Body3
export Body
export Vector3
export unpack
export step!
export GRAVITATIONAL_CONSTANT


import Base.+
import Base.*
import Base.-

const GRAVITATIONAL_CONSTANT = 6.67430e-11

struct Vector3
    x::Float64
    y::Float64
    z::Float64
end
(+)(a::Vector3, b::Vector3) = Vector3(
    a.x + b.x,
    a.y + b.y,
    a.z + b.z
)
(*)(s::Number, u::Vector3) = Vector3(
    s * u.x,
    s * u.y,
    s * u.z
)
(*)(u::Vector3, s::Number) = s * u
(-)(a::Vector3, b::Vector3) = a + -1*b
norm(u::Vector3) = sqrt(u.x^2 + u.y^2 + u.z^2)
unpack(u::Vector3) = [u.x, u.y, u.z]

mutable struct Body 
    r::Vector3
    ṙ::Vector3
    r̈::Vector3
    mass::Float64
    Body(x::Float64, y::Float64, z::Float64) = new(
        Vector3(x,y,z),
        Vector3(0., 0., 0.),
        Vector3(0., 0., 0.),
        1.0
    )
    Body(x::Float64, y::Float64, z::Float64, m::Float64) = new(
        Vector3(x,y,z),
        Vector3(0., 0., 0.),
        Vector3(0., 0., 0.),
        m
    )
    Body(u::Vector3) = new(
        u,
        Vector3(0.,0.,0.),
        Vector3(0.,0.,0.),
        1.0
    )
end

mutable struct Body3
    b1::Body
    b2::Body
    b3::Body
    G::Float64  # Gravitational Constant
    dt::Float64 # Time Differential
    Body3(r1::Body, r2::Body, r3::Body) = new(
        r1, r2, r3,
        GRAVITATIONAL_CONSTANT,
        0.01
    )
    Body3(r1::Body, r2::Body, r3::Body, G::Float64, dt::Float64) = new(
        r1, r2, r3,
        G, dt
    )
end

function body_accel_interaction(b₁::Body, b₂::Body, G::Float64)
    m₁ = b₁.mass;   m₂ = b₂.mass
    r₁ = b₁.r;      r₂ = b₂.r
    (-G * m₂ / norm(r₁ - r₂)^3) * (r₁ - r₂)
end

function advance_acceleration!(prob::Body3)
    b₁ = prob.b1; b₂ = prob.b2; b₃ = prob.b3; G = prob.G
    r̈₁ = body_accel_interaction(b₁, b₂, G) +
         body_accel_interaction(b₁, b₃, G)
    r̈₂ = body_accel_interaction(b₂, b₁, G) +
         body_accel_interaction(b₂, b₃, G)
    r̈₃ = body_accel_interaction(b₃, b₁, G) +
         body_accel_interaction(b₃, b₂, G)

    prob.b1.r̈ = r̈₁
    prob.b2.r̈ = r̈₂
    prob.b3.r̈ = r̈₃

end

function advance_velocity!(prob::Body3)
    dt = prob.dt
    ṙ₁ = prob.b1.ṙ; ṙ₂ = prob.b2.ṙ; ṙ₃ = prob.b3.ṙ
    r̈₁ = prob.b1.r̈; r̈₂ = prob.b2.r̈; r̈₃ = prob.b3.r̈

    ṙ₁ = ṙ₁ + r̈₁ * dt * dt
    ṙ₂ = ṙ₂ + r̈₂ * dt * dt
    ṙ₃ = ṙ₃ + r̈₃ * dt * dt

    prob.b1.ṙ = ṙ₁
    prob.b2.ṙ = ṙ₂
    prob.b3.ṙ = ṙ₃
end

function advance_position!(prob::Body3)
    dt = prob.dt
    r₁ = prob.b1.r; r₂ = prob.b2.r; r₃ = prob.b3.r
    ṙ₁ = prob.b1.ṙ; ṙ₂ = prob.b2.ṙ; ṙ₃ = prob.b3.ṙ

    r₁ = r₁ + ṙ₁ * dt
    r₂ = r₂ + ṙ₂ * dt
    r₃ = r₃ + ṙ₃ * dt

    prob.b1.r = r₁
    prob.b2.r = r₂
    prob.b3.r = r₃
end

function step!(prob::Body3)
    advance_acceleration!(prob)
    advance_velocity!(prob)
    advance_position!(prob)
end

end