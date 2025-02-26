#include <SFML/Graphics.hpp>
#include <variant>
#include "math/Vector2D.hpp"  // Ensure your Vector2D has getX(), getY(), setX(), and setY()

// A simple structure to represent a physics object
struct PhysicsObject {
    Vector2D position;
    Vector2D velocity;
    Vector2D acceleration;
    float radius;

    PhysicsObject(const Vector2D& pos, const Vector2D& vel, const Vector2D& acc, float r)
        : position(pos), velocity(vel), acceleration(acc), radius(r) {}

    // Simple update using semi-implicit Euler integration
    void update(float dt) {
        velocity = velocity + acceleration * dt;
        position = position + velocity * dt;
    }
};

int main() {
    // Create a window using SFML 3.0 VideoMode (using a Vector2u for the dimensions)
    sf::RenderWindow window(sf::VideoMode(sf::Vector2u(800, 600)), "Physics Simulation with SFML");

    // Create a physics object (a ball) with an initial horizontal velocity
    PhysicsObject ball(Vector2D(100.0, 100.0), Vector2D(150.0, 0.0), Vector2D(0.0, 200.0), 20.0f);

    // Create an SFML circle shape to represent the ball
    sf::CircleShape circle(ball.radius);
    circle.setFillColor(sf::Color::Green);
    circle.setOrigin(sf::Vector2f(ball.radius, ball.radius)); // Center the circle

    sf::Clock clock;
    const float restitution = 0.8f; // Bounce coefficient (0.0 = no bounce, 1.0 = perfect bounce)

    while (window.isOpen()) {
        float dt = clock.restart().asSeconds();

        // SFML 3.0 event loop using std::optional and variant-based events
        while (auto event = window.pollEvent()) {
            if (event->is<sf::Event::Closed>())
                window.close();
        }

        // Update the physics for the ball
        ball.update(dt);

        // Get window dimensions for boundary checking
        sf::Vector2u windowSize = window.getSize();

        // Check collision with left and right boundaries
        if (ball.position.getX() - ball.radius < 0) {
            ball.position.setX(ball.radius);
            ball.velocity.setX(-ball.velocity.getX() * restitution);
        }
        if (ball.position.getX() + ball.radius > windowSize.x) {
            ball.position.setX(windowSize.x - ball.radius);
            ball.velocity.setX(-ball.velocity.getX() * restitution);
        }

        // Check collision with top and bottom boundaries
        if (ball.position.getY() - ball.radius < 0) {
            ball.position.setY(ball.radius);
            ball.velocity.setY(-ball.velocity.getY() * restitution);
        }
        if (ball.position.getY() + ball.radius > windowSize.y) {
            ball.position.setY(windowSize.y - ball.radius);
            ball.velocity.setY(-ball.velocity.getY() * restitution);
        }

        // Update the circle's position (convert from your Vector2D to an SFML vector)
        circle.setPosition(sf::Vector2f(static_cast<float>(ball.position.getX()),
                                        static_cast<float>(ball.position.getY())));

        // Clear the window, draw the circle, and display the new frame
        window.clear();
        window.draw(circle);
        window.display();
    }

    return 0;
}

