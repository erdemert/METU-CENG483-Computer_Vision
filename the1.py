from re import A
import numpy as np
from PIL import Image
import sys #First argument will ve the bin number, second is the query number.
#Terminal command : python3 the1.py bin_number choose_query is_histogram_3d grid_size i.e 96
#bin_number = 2,4,8,16,32 etc. 
#choose_query = 1,2, or 3
#is_histogram_3d = 0 or 1
#grid_size = 96 to avoid dividing the image. 48, 24,16 and 12 are possible arguments.


class hist_3_ch() : 
    def __init__(self) :
        self.r_hist = np.zeros(int(sys.argv[1]))
        self.g_hist = np.zeros(int(sys.argv[1]))
        self.b_hist = np.zeros(int(sys.argv[1]))

def l1_normalizer(hist) : 
    out_hist = hist_3_ch()
    out_hist.r_hist = hist.r_hist / np.linalg.norm(hist.r_hist, 1)
    out_hist.g_hist = hist.g_hist / np.linalg.norm(hist.g_hist, 1)
    out_hist.b_hist = hist.b_hist / np.linalg.norm(hist.b_hist, 1)
    return out_hist

def l1_normalizer_3d(hist) : 
    hist = hist / np.linalg.norm(hist, 1)
    return hist

def kulback_div(hist_s, hist_q) : #Takes two histograms and returns the divergence. 
    length = len(hist_s)
    sum = 0
    for i in range (length) : 
        if hist_s[i] == 0: 
            hist_s[i] = 0.000000001
        if hist_q[i] == 0: 
            hist_q[i] = 0.000000001
        sum = hist_q[i] * np.log10(hist_q[i] / hist_s[i]) + sum
    return sum

def calculate_perch_histogram(number_of_bins, img, grid_size) : 
    array_len = int(pow(96/grid_size, 2))
    hist_array = []
    c1 = c2 = 0
    for pieces in range (int(array_len)) : 
        hist_ = hist_3_ch()
        for i in range (grid_size * c1 , grid_size * (c1 + 1)) :   
            for j in range (grid_size * c2,  grid_size * (c2 + 1)) : 
                hist_.r_hist[img[i][j][0] //(256 // int(number_of_bins)) ] += 1
                hist_.g_hist[img[i][j][1] //(256 // int(number_of_bins)) ] += 1
                hist_.b_hist[img[i][j][2] //(256 // int(number_of_bins)) ] += 1
        hist_array.append(l1_normalizer(hist_))
        if c1 == (96//grid_size) - 1 : 
            c1 = 0
            c2 = c2 + 1
        else : 
            c1 = c1 + 1

    return hist_array

def calculate_3ch_histogram(number_of_bins, img, grid_size) : 
    array_len = int(pow(96/grid_size, 2))
    hist_array = []
    c1 = c2 = 0
    for pieces in range (int(array_len)) : 
        bin_number = int(number_of_bins)
        hist_ = np.zeros(pow(bin_number, 3))
        for i in range (grid_size * c1 , grid_size * (c1 + 1)) : 
            for j in range (grid_size * c2, grid_size * (c2 + 1)) : 
                index1 = int(img[i][j][0] // (256 // bin_number)) 
                index2 = int(img[i][j][1] // (256 // bin_number)) 
                index3 = int(img[i][j][2] // (256 // bin_number))
                hist_[index1*bin_number*bin_number + index2 * bin_number + index3 ] = hist_[index1*bin_number*bin_number + index2 * bin_number + index3] + 1
        hist_array.append(l1_normalizer_3d(hist_))
        if c1 == (96//grid_size) - 1 : 
            c1 = 0
            c2 = c2 + 1
        else : 
            c1 = c1 + 1
    return hist_array

np_inst_names = np.loadtxt(
    "InstanceNames.txt", dtype=str)  # Take instance names

all_support_histograms = []

for support_images in np_inst_names : 
    support_image = Image.open("support_96/" + support_images)
    np_img = np.array(support_image)
    if sys.argv[3] == '0' : 
        all_support_histograms.append(calculate_perch_histogram(sys.argv[1], np_img, int(sys.argv[4]) ))
    else : 
        all_support_histograms.append(calculate_3ch_histogram(sys.argv[1], np_img, int(sys.argv[4])))

correct = 0
count = 0
for query_images in np_inst_names :
    count = count + 1
    print("Currently executing query image : " + str(query_images))
    print("Process is at %" + str(count/2))
    query_image = Image.open("query_" + sys.argv[2] +"/" +  query_images)
    np_img_q = np.array(query_image)  # image is converted to numpy array
    if(sys.argv[3] == '0') :
        query_hist = calculate_perch_histogram(sys.argv[1], np_img_q, int(sys.argv[4]))
    else :
        query_hist = calculate_3ch_histogram(sys.argv[1], np_img_q, int(sys.argv[4]))

    div = 9999
    curr_i = 9999

    if sys.argv[3] == '1' : 
        for ith in range(len(all_support_histograms)) : 
            div_f = 0
            for piece in range (pow(96 // int(sys.argv[4]),2)) : 
                div__ = kulback_div(all_support_histograms[ith][piece], query_hist[piece])
                div_f = div_f + div__
            div_f = div_f / pow(96//int(sys.argv[4]), 2 )  
            if div_f < div : 
                curr_i = ith
                div = div__
    else :    
        for ind in range(len(all_support_histograms)) : 
            div_f = 0
            for piece in range ( pow(96 // int(sys.argv[4]),2) ) : 
                div_r = kulback_div(all_support_histograms[ind][piece].r_hist, query_hist[piece].r_hist)
                div_g = kulback_div(all_support_histograms[ind][piece].g_hist, query_hist[piece].g_hist)
                div_b = kulback_div(all_support_histograms[ind][piece].b_hist, query_hist[piece].b_hist)
                div_ = ( div_r + div_g + div_b ) / 3
                div_f = div_f + div_
            div_f = div_f / pow(96//int(sys.argv[4]), 2 )  
            if div_f < div : 
                curr_i = ind
                div = div_f
    
    if query_images == np_inst_names[curr_i] : 
        correct = correct + 1

print("Top 1 accuracy is: " + str(correct/len(np_inst_names)) + "! " +str(sys.argv[1]) + " bins are used. " + str(sys.argv[4]) + " grid is used at query " +
str(sys.argv[2]))

