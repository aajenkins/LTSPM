# @Author: Jenkins Alec <alec>
# @Date:   2017-01-28T16:36:20-08:00
# @Project: LTSPM analysis
# @Last modified by:   alec
# @Last modified time: 2017-01-28T16:42:45-08:00

# find position of contours in an array and return a list of contours,
# each as an ordered list of coordinates. takes binary image array and
# starting positions of each contour

def find_cotours(image, starting_points):
