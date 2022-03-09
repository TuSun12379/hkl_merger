import hkl_IO as hkl_io
import numpy as np
import matplotlib.pyplot as plt

def plot_theta_intensity(theta_list, inte_list, show_hkl=False):
    """To plot intensity against 2theta.
    """
    theta_list, inte_list = np.array(theta_list), np.array(inte_list).reshape(-1, 5)
    plt.figure(figsize=(22, 8), dpi=100)
    plt.xlabel("2 theta, ${\lambda}=0.02508\,{\AA}$")
    plt.ylabel("Intensity")
    plt.scatter(theta_list, inte_list[:, 3], s = 30, c='green', marker='o', alpha=0.7)
    if show_hkl:
        for i in range(len(theta_list)):
            plt.text(theta_list[i], inte_list[i, 3],
                    "{:.0f} {:.0f} {:.0f}".format(inte_list[i, 0],
                    inte_list[i, 1], inte_list[i, 2]))
    plt.show()

def plot_theta_intensity2(theta_list1, inte_list1, theta_list2, inte_list2, show_hkl=False):
    """To plot intensity against 2theta.
    """
    theta_list1, inte_list1 = np.array(theta_list1), np.array(inte_list1).reshape(-1, 5)
    theta_list2, inte_list2 = np.array(theta_list2), np.array(inte_list2).reshape(-1, 5)
    plt.figure(figsize=(22, 8), dpi=100)
    axes = plt.subplot(111)
    label_nonabset = axes.scatter(theta_list1, inte_list1[:, 3], s=30, c='green', marker='o', alpha=0.7)
    label_absent = axes.scatter(theta_list2, inte_list2[:, 3], s=30, c='red', marker='o', alpha=0.7)
    plt.xlabel('2 theta, ${\lambda}=0.02508\,{\AA}$')
    plt.ylabel('Intensity')
    axes.legend((label_nonabset, label_absent),
                ('hkl-1',
                'hkl-2'), loc=1)
    if show_hkl:
        for i in range(len(theta_list1)):
            plt.text(theta_list1[i], inte_list1[i, 3],
                    "{:.0f} {:.0f} {:.0f}".format(inte_list1[i, 0], inte_list1[i, 1], inte_list1[i, 2]))
        for j in range(len(theta_list2)):
            plt.text(theta_list2[j], inte_list2[j, 3],
                    "{:.0f} {:.0f} {:.0f}".format(inte_list2[j, 0], inte_list2[j, 1], inte_list2[j, 2]))
    plt.show()

if __name__ == '__main__':
    pass
