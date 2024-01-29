import yt

ds = yt.load("../Data_000150")

slc = yt.SlicePlot(
    ds,
    "z",
    [("gas", "density"), ("gas", "kT"), ("gas", "pressure"), ("gas", "velocity_y")],
)
slc.set_unit(("gas", "velocity_y"), "c")
slc.save()
