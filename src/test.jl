f(x,y) = (3x + y^2) * abs(sin(x) + cos(y))

x = 0:0.01:2
y = 0:0.01:1
z = [f(i,j) for i in x, j in y]' # この転置を忘れるとデータが矛盾し、グラフが変になるので要注意．
print(size(x))
print(size(y))
print(size(z))
plot(x,y,z, st=:surface)