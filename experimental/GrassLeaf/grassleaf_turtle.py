import turtle as t

t.shape("turtle")

x = 1
t.speed(10)
colors = ["red", "orange", "yellow", "green", "blue", "purple"]

for i in range(100):
    for i in colors:
        t.forward(x * 0.3)
        t.left(60)
        t.color(i)
        t.right(30.5)
        x = x + 1
