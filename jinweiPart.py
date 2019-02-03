import sys
import numpy as np
from PNM import *

def CreateAndSavePFM(out_path):
    width = 512
    height = 512
    numComponents = 3

    img_out = np.empty(shape=(width, height, numComponents), dtype=np.float32)

    for y in range(height):
        for x in range(width):
            img_out[y,x,:] = 1.0

    writePFM(out_path, img_out)

def LoadAndSavePPM(in_path, out_path):
    img_in = loadPPM(in_path)
    img_out = np.empty(shape=img_in.shape, dtype=img_in.dtype)
    height,width,_ = img_in.shape # Retrieve height and width
    for y in range(height):
        for x in range(width):
            img_out[y,x,:] = img_in[y,x,:] # Copy pixels

    writePPM(out_path, img_out)

def LoadAndSavePFM(in_path, out_path):
    img_in = loadPFM(in_path)
    img_out = np.empty(shape=img_in.shape, dtype=img_in.dtype)
    height,width,_ = img_in.shape # Retrieve height and width
    for y in range(height):
        for x in range(width):
            img_out[y,x,:] = img_in[y,x,:] # Copy pixels

    writePFM(out_path, img_out)

def LoadPPMAndSavePFM(in_path, out_path):
    img_in = loadPPM(in_path)
    img_out = np.empty(shape=img_in.shape, dtype=np.float32)
    height,width,_ = img_in.shape
    for y in range(height):
        for x in range(width):
            img_out[y,x,:] = img_in[y,x,:]/255.0

    writePFM(out_path, img_out)

def LoadPFMAndSavePPM(in_path, out_path):
    img_in = loadPFM(in_path)
    img_out = np.empty(shape=img_in.shape, dtype=np.float32)
    height,width,_ = img_in.shape
    for y in range(height):
        for x in range(width):
            img_out[y,x,:] = img_in[y,x,:] * 255.0

    writePPM(out_path, img_out.astype(np.uint8))


def MapingLatlong(input_path, out_path):
    """
    :param input_path:
    :param out_path:
    :return:
    """
    latlong = loadPFM(input_path)
    print('shape of the latlong', latlong.shape)

    y_length, x_length, channel = latlong.shape

    # initialize the output image
    diameter = 511
    radius = float(diameter)/2
    v = [0.0, 0.0, 1.0] #viewing vector
    light_Probe = np.empty(shape=(diameter, diameter, channel), dtype=np.float32)
    print('shape of the light Prob', light_Probe.shape)

    # do the mapping for each pixels
    for h in range(0, diameter): # height of the image
        for w in range(0, diameter): # weight of the image
            dis = np.sqrt(np.power(h-radius, 2) + np.power(w-radius, 2))
            # distance between each pixel to the center
            if dis >= radius:
                for ch in range(3):
                    light_Probe[h][w][ch] = 0.0
            else:
                #map to polar coordinate system
                n = cartToPolar(h,w,radius)
                # this n has been normalised, reflection vector
                r = (2*(np.dot(n, v))*n - v)
                # the Cartesian coordinates after conversion
                a, b = polarToCart(r, y_length, x_length)
                # mapping back the sphere
                for ch in range(channel):
                    if latlong[b][a][ch] > 1.0:
                        latlong[b][a][ch] = 1.0
                    light_Probe[h][w][ch] = latlong[b][a][ch]

    writePFM(out_path, light_Probe)

def cartToPolar(height, weight, radius):
    """
    :param height:
    :param weight:
    :param radius:
    :return:
    """
    y = radius - height
    x = weight - radius
    z = np.sqrt(np.power(radius, 2) - np.power(y, 2) - np.power(x, 2))

    if z == 0:
        norm = [x / radius, y / radius, 0]
    else:
        phi = np.arctan2(y, x)  # phi
        theta = np.arccos(z / radius)  # theta
        norm = np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])

    return norm
def polarToCart(rVector, y_length, x_length):
    """
    :param rVector:
    :param y_length:
    :param x_length:
    :return:
    """
    # chang the angle phi from (-pi , pi) to (0, 2pi)
    x, y, z = rVector[0], rVector[1], rVector[2]
    R_y = np.arctan2(x, z) + np.pi
    R_projectZ = np.arccos(y)

    a = int(((R_y) / (2 * np.pi)) * x_length)
    b = int((R_projectZ / np.pi) * y_length)

    return a, b

def GammaCorrection(img_path, gamma, stop, out_path):
    """
    :param img_path:
    :param stop:
    :param gamma:
    :param out_path:
    :return:
    """
    img = loadPFM(img_path)
    height, width, channel = img.shape

    G = np.empty(shape=img.shape, dtype=img.dtype)


    for h in range(height):
        for w in range(width):
            for ch in range(channel):
                img[h][w][ch] *= stop # scaling factorR
                img[h][w][ch] = np.power(img[h][w][ch], 1/gamma)  # gamma correction
                if img[h][w][ch] > 1.0:
                    G[h][w][ch] = 255
                else:
                    G[h][w][ch] = img[h][w][ch]*255
    writePPM(out_path, G.astype(np.uint8))

if '__main__' == __name__:
    # LoadAndSavePFM('grace_latlong.pfm', 'test.pfm')
    # LoadAndSavePPM('9.ppm', 'test.ppm')
    # LoadPFMAndSavePPM('test.pfm', 'grace.ppm')
    # LoadPPMAndSavePFM('test.ppm', '9.pfm')
    # MapingLatlong('../UrbanProbe/urbanEM_latlong.pfm','../UrbanProbe/result2.pfm')
    gamma = [1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8]
    stop = [3,5,7,9]
    for g in [2.5]:
        for s in stop:
            GammaCorrection('../UrbanProbe/result2.pfm', g, s, '../UrbanProbe/result2AfterGamma'+str(g)+'Stop'+str(s)+'.pfm' )
