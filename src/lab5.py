from PIL import Image
import numpy as np

n = 300

input_array = np.ones((n,n,3)) * 255

model_image = Image.new('RGB', (n, n), color = 'red')
model_image.save('../Images/model_image.png')

input_image = Image.fromarray(input_array.astype('uint8')).convert('RGB')
input_image.save('../Images/input_image.png')


def getRandomDirection():
    r1 = np.random.uniform()
    r2 = np.random.uniform()
    x = np.cos(2*np.pi*r1)*np.power(r2*(1-r2), 0.5)
    y = np.sin(2*np.pi*r1)*np.power(r2*(1-r2), 0.5)
    z = 1 - 2*r2
    return np.array([x, y, z])
def colorMatching(model_image, input_image, num_iter=1000):
    model = Image.open('../Images/model_image.png')
    input = Image.open('../Images/input_image.png')


    for i in range(num_iter):
        print("Iteration: ", i + 1)
        v = getRandomDirection()
        projI = []
        projM = []
        for j in range(n*n):
            x = j%n
            y = j//n
            rm, gm, bm = model.getpixel((x, y))[0], model.getpixel((x, y))[1], model.getpixel((x, y))[2]
            ri, gi, bi = input.getpixel((x, y))[0], input.getpixel((x, y))[1], input.getpixel((x, y))[2]
            pM = np.dot(np.array([rm, gm, bm]), v)
            pI = np.dot(np.array([ri, gi, bi]), v)
            projI.append((pI, j))
            projM.append((pM, j))

        projI.sort(key=lambda x: x[0])
        projM.sort(key=lambda x: x[0])

        for j in range(n*n):
            dot_Input, i_input = projI[j]
            dot_Model, i_model = projM[j]
            x = i_input%n
            y = i_input//n
            old_value = input.getpixel((x, y))
            new_value = old_value + (dot_Model - dot_Input) * v
            new_value = new_value.astype('uint8')
            new_value = np.clip(new_value, 0, 255)
            input.putpixel((x, y), (new_value[0], new_value[1], new_value[2]))

    input.save("input_image_transformed.png")

colorMatching(model_image, input_image)
