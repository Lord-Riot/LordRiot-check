CRT-RSA体系的攻击，当p，q不均衡（p < q），且dq较小时，可以通过多元CopperSmith攻击，这里给的bound非常宽松，只是卡了最初的二元格，使用Bleichenbacher-May attack即可恢复p，最新的paper其实是可以做到p取0.5的，但m太大，考虑到复杂度问题，于是就用了这个0.35的。