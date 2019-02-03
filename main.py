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


######################################### CODE START FROM HERE #########################################


# center weight function
def Gaussian_weight(x):
    mu = .5
    sig = .1
    return np.exp(-(np.power((x-mu), 2.)) / (2.*np.power(sig, 2.)))


# create assemble HDR
def Assemble_HDR():

    PFM_NUM = 7
    # load pfm into Z
    Z = []
    # create relative exposure
    delta_t = []

    for i in range(PFM_NUM):
        img_in = loadPFM('../Office/Office{}.pfm'.format(i+1))
        Z.append(img_in)
        delta_t.append(2.**(i*2.))
    height, width, channel = Z[0].shape

    # cal intensity
    Z_intensity = []
    F = np.empty(shape=(height, width, channel))

    for i in range(PFM_NUM):
        this_intensity = np.empty(shape=(height, width))
        for h in range(height):
            for w in range(width):
                # transfer RGB to intensity
                this_intensity[h, w] = Z[i][h, w, 0]*.3 + Z[i][h, w, 1]*.6 + Z[i][h, w, 2]*.1

        Z_intensity.append(this_intensity)


    # for each pixel and channel calculate weighted avg
    for h in range(height):
        for w in range(width):
            for ch in range(channel):
                F[h, w, ch] = 0.
                total_weight = 0.

                for i in range(PFM_NUM):
                    # good pixel
                    if (Z_intensity[i][h, w] >= 0.005) and (Z_intensity[i][h, w] <= .92):
                        this_weight = Gaussian_weight(Z_intensity[i][h, w])
                        F[h, w, ch] += np.log((1./delta_t[i]) * Z[i][h, w, ch]) * this_weight
                    # bad pixel
                    else:
                        this_weight = 0
                    total_weight += this_weight

                F[h, w, ch] = np.exp(F[h, w, ch]/total_weight)


    # calculate brightest and dimmest point
    upper_f = 0.
    lower_f = 1.
    for h in range(height):
        for w in range(width):
                intensity = F[h, w, 0]*.3 + F[h, w, 1]*.6 + F[h, w, 2]*.1
                if intensity  < lower_f:
                    lower_f = intensity
                if intensity > upper_f:
                    upper_f = intensity

    #calculate radius
    rad = upper_f/lower_f

    print(upper_f, lower_f, rad)

    writePFM('../AssemblePFM.pfm', F)


def linear_scale(F):

    height, width, channel = F.shape

    upper_f = 0.
    for h in range(height):
        for w in range(width):
            intensity = F[h, w, 0] * .3 + F[h, w, 1] * .6 + F[h, w, 2] * .1
            if intensity > upper_f:
                upper_f = intensity
    print upper_f

    for h in range(height):
        for w in range(width):
            for ch in range(3):
                F[h, w, ch] = F[h, w, ch]/upper_f
                if F[h, w, ch] > 1:
                    F[h, w, ch] = 1

    return F


def exp_scale(F, stop):

    height, width, channel = F.shape

    for h in range(height):
        for w in range(width):
            for ch in range(3):
                for t in range(stop):
                    F[h, w, ch] *= 2
                if F[h, w, ch] > 1:
                    F[h, w, ch] = 1

    return F


def gamma_scale(F, gam, stop):

    height, width, channel = F.shape
    for h in range(height):
        for w in range(width):
            for ch in range(channel):
                for i in range(stop):
                    F[h, w, ch] = F[h, w, ch] ** (1. / gam)


    return F


if '__main__' == __name__:

    # Assemble_HDR()

    P = loadPFM('../AssemblePFM.pfm')
    F = P[:]
    height, width, channel = F.shape
    # print F

    gam = 2.6
    exp_stop = 1
    gam_stop = 2

    F = linear_scale(F)
    F = exp_scale(F, exp_stop)
    F = gamma_scale(F, gam, gam_stop)

    for h in range(height):
        for w in range(width):
            for ch in range(3):
                F[h, w, ch] = np.ceil(F[h, w, ch] * 255)
                if F[h, w, ch] > 255:
                    F[h, w, ch] = 255

    # writePFM('../HDR_expstop_{}_gamma{}_gamstop{}.pfm'.format(exp_stop, gam, gam_stop), F)
    writePPM('../new/HDR_expstop_{}_gamma{}_gamstop{}.ppm'.format(exp_stop, gam, gam_stop), F.astype(np.uint8))


