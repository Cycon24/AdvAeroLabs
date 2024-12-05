from PIL import Image
import numpy as np
import matplotlib.pyplot as plt

# set up arrays of the velocities and degrees studied
vs = np.array([150., 200., 250.])
degs = np.array([0., 1.5, 3., 5., 7., 9., 11.])

# create dict for storing data
image_data = {}
for v in vs:
    for deg in degs:
        deg_str = f"{deg:.0f}" if deg.is_integer() else f"{deg:.1f}"
        file_path = rf"Lab_05\ONERA M6 - Surface Pressure\V{v:.0f}_{deg_str}deg.tif"
        
        # open image
        image = Image.open(file_path)
        image_data[(v, deg)] = np.array(image)

# wing pixel locs
FO = [85.5, 112.6]
FI = [778.1, 35.5]
AO = [19.3, 314.9]
AI = [783.2, 341.1]

# dims of onera M6
cr = 5.39
ct = 3.03
lam = 0.562

# measure distance between pixel locations
ctpix = np.sqrt((FO[0] - AO[0])**2 + (FO[1] - AO[1])**2)
# crpix = np.sqrt((FI[0] - AI[0])**2 + (FI[1] - AI[1])**2) ignoring this measurement because cr is measured in fuselage
crpix = ctpix/lam                                        # using taper ratio to estimate instead for better agreement

# calculate pixels/inch and average to get scale of data
ctrat = ctpix/ct
crrat = crpix/cr
ratio = np.mean([ctrat, crrat])
print(f"Image Scale: {ratio:.2f} pix/in")

# convert the image to array
image_array = image_data[(150, 0)]
height, width = image_array.shape

# scale axes
x_axis = np.arange(width)/ratio
y_axis = np.arange(height)/ratio


# calculate Cp
def Cp(V, P):
    # static pressure taken from data at 0 aoa
    Ps = 1.0059E5
    
    # density assumption
    rho = 1.225
    
    # dynamic pressure
    V /= 3.281
    q = 0.5*rho*V**2
    
    Cp = (P - Ps)/q
    return Cp


# cp calc function
def cpsC(V, image_array):
    cps = np.empty_like(image_array)
    
    # loop through vals and calc Vp with function
    for i in range(len(image_array[:, 0])):
        for k in range(len(image_array[0, :])):
            val = image_array[i, k]
            if val == val:
                cps[i, k] = Cp(V, val)

    return cps


# plotting function
def contourPlot(V, deg):
    image_array = image_data[(V, deg)]
    cps = cpsC(V, image_array)
    plt.figure(f"V={V:.0f} ft/s @ {deg:.1f} deg")
    levels = np.linspace(-6., 3.5, 100)
    contour = plt.contourf(x_axis, y_axis, cps, levels=levels, cmap='magma')
    colorbarticks = np.arange(-6, 3.51, 0.5)
    plt.colorbar(contour, ticks=colorbarticks, label="Pressure Coefficient $C_P$")
    plt.xlabel("Width (in.)"); plt.ylabel("Height (in.)")
    plt.title(f"V={V:.0f} ft/s @ {deg:.1f} deg")
    plt.gca().invert_yaxis()
    plt.gca().set_aspect("equal")
    plt.ylim(6, 0); plt.xlim(0, 12)
    plt.grid(); plt.tight_layout()
    deg_str = f"{deg:.0f}" if deg.is_integer() else f"{deg:.1f}"
    plt.savefig(rf"Lab_05/figs/V{V:.0f}_{deg_str}deg.png", dpi=300)


# loop through and plot everything
for V in vs:
    for deg in degs:
        contourPlot(V, deg)
        
# pull out pressure vals at specific locs
locs = np.array([0.2, 0.44, 0.6, 0.8, 0.9])

# select points along midpoint of wing
root = np.array([11., 3.])*ratio      # root location (pix)
tip = np.array([0.5, 3.])*ratio       # tip loc (pix)
b = root[0] - tip[0]                  # span (pix)

# get analysis locations based on span
locsx = np.array([b*(1 - loc) for loc in locs])
locsy = np.array([root[1] for _ in range(len(locs))])

# pull out pressure vals at location
def pressPlot(V, deg, locsx, locsy):
    # get cp data
    image_array = image_data[(V, deg)]
    cps = cpsC(V, image_array)
    
    # get press at each location
    press = np.array([np.mean(cps[int(locsy[i]) - 5:int(locsy[i] + 5), int(locsx[i]) - 5: int(locsx[i] + 5)]) for i in range(len(locsx))])
    plt.plot(locs, press, label=f"{deg:.1f}$^\circ$", linewidth=2)

# plot contours
plt.figure()
[pressPlot(150, deg, locsx, locsy) for deg in degs]
plt.legend(loc="upper left", fontsize=12)
plt.xlim(0, 1); plt.ylim(-4, 0)
plt.xlabel("Spanwise Location (y/b)"); plt.ylabel("Pressure Coefficient $C_P$")
plt.title("V=150 ft/s")
plt.gca().invert_yaxis()
plt.grid(); plt.tight_layout()
plt.savefig(r"Lab_05/figs/V150prescont.png", dpi=300)

plt.show()