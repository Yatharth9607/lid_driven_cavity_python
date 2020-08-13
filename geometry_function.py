# geometry function


def geometry_function(nx, ny):
    l = 0.5
    h = 0.5

    ipcv = nx
    jpcv = ny

    dx = l / ipcv
    dy = h / jpcv

    return ipcv, jpcv, dx, dy, l, h
