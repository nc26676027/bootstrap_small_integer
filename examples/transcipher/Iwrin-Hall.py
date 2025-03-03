import numpy as np
from scipy.stats import irwinhall

# 设置参数 n
n = 2**4 + 1
K = 7
N = 16384
# 设置 loc 和 scale 参数
loc = -0.5 * n
scale = 1 # [-0.5,0.5] interval length 1

cdf_K = irwinhall.cdf(K, n, loc=loc, scale=scale)
print("cdf_K:", cdf_K)

# Check if cdf_K is sufficiently less than 1 for valid computation
p = max(0, 2*cdf_K-1)  # Ensure non-negative and less than or equal to 1

if p < 1:
    res = 1 - p**N
    probability_log2 = np.log2(res)
else:
    probability_log2 = float('-inf')  # Log(0) situation

print("Pr(>K):", probability_log2)
