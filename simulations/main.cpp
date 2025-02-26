#include <SFML/Graphics.hpp>
#include "math/Vector2D.hpp"  // Adjust path if necessary

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
    // Use SFML 3.0 VideoMode that takes sf::Vector2u
    sf::RenderWindow window(sf::VideoMode(sf::Vector2u(800, 600)), "Physics Simulation with SFML");

    // Create a physics object (e.g., a falling ball)
    PhysicsObject ball(Vector2D(100.0, 100.0), Vector2D(0.0, 0.0), Vector2D(0.0, 200.0), 20.0f);

    // Create an SFML shape to represent the ball
    sf::CircleShape circle(ball.radius);
    circle.setFillColor(sf::Color::Green);
    circle.setOrigin(sf::Vector2f(ball.radius, ball.radius)); // Center the circle

    sf::Clock clock;

    while (window.isOpen()) {
        // Calculate elapsed time
        float dt = clock.restart().asSeconds();

        // Updated event loop using std::optional
        while (auto event = window.pollEvent()) {
            if (event->is<sf::Event::Closed>())
                window.close();
        }

        // Update the physics object
        ball.update(dt);

        // Update the circle's position using getters from Vector2D
        circle.setPosition(sf::Vector2f(static_cast<float>(ball.position.getX()),
                                        static_cast<float>(ball.position.getY())));

        // Clear and draw
        window.clear();
        window.draw(circle);
        window.display();
    }

    return 0;
}

