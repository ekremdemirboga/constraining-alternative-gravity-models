# importing the modules needed
import cv2
import numpy as np
  
# Reading the image
image = cv2.imread('MRprob_M28.jpeg')
  
  
kernel2 = np.ones((10,10),np.float32)/100
Z = cv2.filter2D(src=image, ddepth=-1, kernel=kernel2)
# Shoeing the original and output image


  
cv2.waitKey()
cv2.destroyAllWindows()