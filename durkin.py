import customtkinter as ctk
from tkinter import *
import matplotlib.pyplot
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,  
NavigationToolbar2Tk)
from pycraf import pathprof, conversions
from astropy import units as u
import numpy as np
from numpy.linalg import norm
from scipy.constants import speed_of_light as speed

ctk.set_default_color_theme("dark-blue")


class Durkin(ctk.CTkToplevel):
    
    APP_NAME = "Result"
    WIDTH = 1000
    HEIGHT = 900

    def __init__(self, tx, rx, f, gt, gr, htx, hrx):
        super().__init__()
        self.tx = tx
        self.rx = rx
        self.f = f
        self.gt = gt
        self.gr = gr
        self.htx = htx
        self.hrx = hrx

        self.title(Durkin.APP_NAME)
        self.geometry(str(Durkin.WIDTH) + "x" + str(Durkin.HEIGHT))
        self.minsize(Durkin.WIDTH, Durkin.HEIGHT)

        loss_label= ctk.CTkLabel(master=self, text=str(), font=("Helvetica Bold",18))
        loss_itu_label = ctk.CTkLabel(master=self, text=str(), font=("Helvetica Bold",18))
        loss_label.pack(pady=(10,0))
        loss_itu_label.pack()
        self.tabview = ctk.CTkTabview(master=self, corner_radius=10, 
                         width=850, height=800)
        self.tabview.add("Elevation")
        self.tabview.add("Terrain")
        self.tabview.add("Propagation")
        self.tabview.pack()

        self.fig_elevation = Figure(figsize=(10, 10), dpi=100)
        self.fig_terrain = Figure(figsize=(10, 10), dpi=100)
        self.fig_propagation =  Figure(figsize=(10, 10), dpi=100)
        self.elevation_plot = self.fig_elevation.add_subplot(111)
        self.terrain_plot = self.fig_terrain.add_subplot(111)
        self.propagation_plot = self.fig_propagation.add_subplot(111)

            # allow download of missing SRTM data:
        pathprof.SrtmConf.set(download='missing', srtm_dir="./srtm")

        freq = f/1000 * u.GHz
        f = f*(10**6)
        lon_tx, lat_tx = tx[0]* u.deg,  tx[1] * u.deg
        lon_rx, lat_rx = rx[0]* u.deg, rx[1] * u.deg
        hprof_step = 100 * u.m
        omega = 0. * u.percent  # fraction of path over sea
        temperature = 290. * u.K
        pressure = 1013. * u.hPa
        time_percent = 2 * u.percent  # see P.452 for explanation
        h_tx, h_rx = htx * u.m, hrx * u.m
        G_t, G_r = gt * conversions.dBi, gr * conversions.dBi
        zone_t, zone_r = pathprof.CLUTTER.UNKNOWN, pathprof.CLUTTER.UNKNOWN
        map_size_lon, map_size_lat = 1.5 * u.deg, 1.5 * u.deg
        map_resolution = 3. * u.arcsec
        
        lons, lats, heightmap = pathprof.srtm_height_map(
        lon_tx, lat_tx,
        map_size_lon, map_size_lat,
        map_resolution=map_resolution)

        terrainmap = heightmap.to(u.m).value
        _lons = lons.to(u.deg).value
        _lats = lats.to(u.deg).value
        vmin, vmax = -20, 170 #limits
        terrain_cmap, terrain_norm = pathprof.terrain_cmap_factory(vmin=vmin, vmax=vmax)
        terrainmap[terrainmap < 0] = 0.51  # fix for coastal region
        ax = self.fig_terrain.add_axes((0., 0., 1.0, 1.0))
        cbax = self.fig_terrain.add_axes((0., 0., 1.0, .02))
        cim = ax.imshow(
            terrainmap,
            origin='lower', interpolation='nearest',
            cmap=terrain_cmap, norm=terrain_norm,
            extent=(_lons[0], _lons[-1], _lats[0], _lats[-1]))
        cbar = self.fig_terrain.colorbar(cim, cax=cbax, orientation='horizontal')
        
        ax.set_aspect(abs(_lons[-1] - _lons[0]) / abs(_lats[-1] - _lats[0]))

        ctics = np.arange(0, vmax, 50)
        cbar.set_ticks(ctics)
        cbar.ax.set_xticklabels(map('{:.0f} m'.format, ctics), color='k')
        cbar.set_label(r'Height (amsl)', color='k')
        cbax.xaxis.tick_top()
        cbax.xaxis.set_label_position('top')

        ax.set_xlabel('Longitude [deg]')
        ax.set_ylabel('Latitude [deg]')

        map_size_lon, map_size_lat = 0.5 * u.deg, 0.5 * u.deg
        map_resolution = 10. * u.arcsec

        hprof_cache = pathprof.height_map_data(
        lon_tx, lat_tx,
        map_size_lon, map_size_lat,
        map_resolution=map_resolution,
        zone_t=zone_t, zone_r=zone_r,
        )

        results = pathprof.atten_map_fast(
        freq,
        temperature,
        pressure,
        h_tx, h_rx,
        time_percent,
        hprof_cache,
        )
        lons = hprof_cache['xcoords']
        lats = hprof_cache['ycoords']
        total_atten = results['L_b'].value

        vmin, vmax = -5, 195
        ax = self.fig_propagation.add_axes((0., 0., 1.0, 1.0))
        cbax = self.fig_propagation.add_axes((0., 0., 1.0, .02))
        cim = ax.imshow(
            total_atten,
            origin='lower', interpolation='nearest', cmap='inferno_r',
            vmin=vmin, vmax=vmax,
            extent=(lons[0], lons[-1], lats[0], lats[-1]),
            )
        cbar = self.fig_propagation.colorbar(
            cim, cax=cbax, orientation='horizontal'
            )
        ax.set_aspect(abs(lons[-1] - lons[0]) / abs(lats[-1] - lats[0]))
        ctics = np.arange(0, vmax, 30)
        cbar.set_ticks(ctics)
        cbar.ax.set_xticklabels(map('{:.0f} dB'.format, ctics), color='w')
        cbar.set_label(r'Path propagation loss', color='w')
        cbax.xaxis.tick_top()
        cbax.tick_params(axis='x', colors='w')
        cbax.xaxis.set_label_position('top')

        ax.set_xlabel('Longitude [deg]')
        ax.set_ylabel('Latitude [deg]')

        (lons,lats,distance,distances,heights,bearing,back_bearing,back_bearings) \
        = pathprof.srtm_height_profile(lon_tx, lat_tx, lon_rx, lat_rx, hprof_step)
        pprop = pathprof.PathProp(  
            freq,
            temperature, pressure,
            lon_tx, lat_tx,
            lon_rx, lat_rx,
            h_tx, h_rx,
            hprof_step,
            time_percent,
            zone_t=zone_t, zone_r=zone_r,
            )
        
        self.loss_itu = pathprof.loss_complete(pprop, G_t, G_r)[-1].value
        self.distance = distance.to(u.m).value

        distances = distances.to(u.km).value
        heights = heights.to(u.m).value
        tx = [distances[0], htx+heights[0]]
        rx = [distances[-1], hrx+heights[-1]]

        self.edges, self.loss = self.durkin(tx, rx, distances, heights)
        loss_label.configure(text=f"Path Loss: {round(self.loss,3)} dB")
        loss_itu_label.configure(text=f"Path Loss (ITU P.452): {round(self.loss_itu,3)} dB")

        elevation = FigureCanvasTkAgg(self.fig_elevation, 
                               master = self.tabview.tab("Elevation"))
        elevation.draw()
        toolbar = NavigationToolbar2Tk(elevation, self.tabview.tab("Elevation")) 
        toolbar.update() 
        elevation.get_tk_widget().pack()

        terrain = FigureCanvasTkAgg(self.fig_terrain, 
                               master = self.tabview.tab("Terrain"))
        terrain.draw()
        terrain.get_tk_widget().pack()

        propagation = FigureCanvasTkAgg(self.fig_propagation, 
                               master = self.tabview.tab("Propagation"))
        propagation.draw()
        propagation.get_tk_widget().pack()


        self.elevation_plot.plot(distances, self.line_of_sight(tx,rx)(distances), "k", linestyle="dotted")
        self.elevation_plot.plot(tx[0], tx[1]+2, 'bv')
        self.elevation_plot.plot(rx[0], rx[1]+2, 'gv')
        self.elevation_plot.annotate("Tx", xy=(tx[0], tx[1]+5), color="b")
        self.elevation_plot.annotate("Rx", xy=(rx[0], rx[1]+5), color="g")
        self.elevation_plot.plot(distances, heights, 'k-')
        self.elevation_plot.set_xlabel('Distância [km]')
        self.elevation_plot.set_ylabel('Altura [m]')


    def get_edge_coordinates(self, tx, rx, x, y):
        los = self.line_of_sight(tx, rx)
        edges_interval = np.where(y > los(x))[0]
        angles = []
        if len(edges_interval) > 1:
            if x[edges_interval[0]] == tx[0]:
                edges_interval = edges_interval[1:]
            for coord in edges_interval:
                txrx = np.subtract(rx,tx)
                txedge = np.subtract([x[coord], y[coord]],tx)
                angle = np.dot(txrx,txedge)/(norm(txrx)*norm(txedge))
                angles.append(np.arccos(angle))
            coordinate = edges_interval[np.argmax(angles)]
            return [x[coordinate], y[coordinate]]
        return None

    def line_of_sight(self, tx, rx):
        los = lambda x: (rx[1]-tx[1])/(rx[0]-tx[0])*x + (rx[0]*tx[1]-rx[1]*tx[0])/(rx[0]-tx[0])
        return los

    def calculate_edge_points(self,tx, rx, x, y):
        if tx:
            edge = self.calculate_edge_points(self.get_edge_coordinates(tx, rx, x, y), rx, x, y)
            if edge:
                return tx+edge
            return tx

    def get_edge_points(self, tx, rx, x, y):
        edges = self.calculate_edge_points(tx, rx, x, y)
        edges = np.array(edges[2:]).reshape(-1,2).tolist()
        bullington = None
        for edge in edges:
            self.elevation_plot.plot(*edge, "ro")
        if len(edges) > 3:
            bullington = self.bullington(tx, rx, edges[1:-1], x)
            edges[1:-1] = [bullington]

        for i in edges:
            los = self.line_of_sight(i,rx)(x)
            if i == bullington:
                xaxis = np.argmin(np.abs(x-bullington[0]))
            else:
                xaxis = np.argmin(np.abs(y-los))
            self.elevation_plot.plot(x[xaxis:], los[xaxis:], 'k', linestyle="dotted")
        return edges

    def fresnel_diffraction_parameter(self, tx, edge, rx):
        tx, edge, rx = map(np.array, [tx, edge, rx])
        h = np.cross(tx-edge, rx-tx)/norm(rx-tx)
        d1 = np.abs(norm(tx-edge)*np.dot(rx-tx,edge-tx)/(norm(rx-tx)*norm(edge-tx)))
        d2 = np.abs(norm(edge-rx)*np.dot(tx-rx,tx-edge)/(norm(tx-rx)*norm(tx-edge)))
        return h*np.sqrt((2*(self.f)*(d1+d2))/(speed*d1*d2))
    

    def diffraction_gain(self, points):
        total_diffraction = 0
        for slide in range(len(points)-2):
            tx, edge, rx = points[slide:slide+3]
            v = self.fresnel_diffraction_parameter(tx, edge, rx)
            if v <= -0.806:
                continue
            elif v >= -0.806 and v <= 0:
                total_diffraction += 20*np.log10(0.5-0.62*v)
            elif v >= 0 and v <= 1:
                total_diffraction += 20*np.log10(0.5*np.exp(-0.95*v))
            elif v >= 1 and v <= 2.4:
                total_diffraction += 20*np.log10(0.4-np.sqrt(0.1184-(0.38-0.1*v)**2))
            else:
                total_diffraction += 20*np.log10(0.225/v)
        return total_diffraction

    def loss(self, tx, rx, edges):
        freespace = -10*np.log10(self.gt*self.gr*((speed/self.f)/(4*np.pi*self.distance))**2)
        flatearth = 40*np.log10(self.distance) - \
        (10*np.log10(self.gt+10*np.log10(self.gr)+20*np.log10(self.htx)+20*np.log10(self.hrx)))
        diffraction_loss = self.diffraction_gain(edges)
        return max(freespace, flatearth)-diffraction_loss


    def bullington(self, tx, rx, edges, distances):
        angles_rx = []
        angles_tx = []
        for edge in edges:
            angles_tx.append(np.arctan2(np.abs(edge[1]-tx[1]), np.abs(edge[0]-tx[0]))+2*np.pi)
            angles_rx.append(np.arctan2(np.abs(edge[1]-rx[1]), np.abs(edge[0]-rx[0]))+2*np.pi)
        edge_tx_coordinate = edges[np.argmax(angles_tx)]
        edge_rx_coordinate = edges[np.argmax(angles_rx)]

        los_edge_tx = self.line_of_sight(tx, edge_tx_coordinate)(distances)
        los_edge_rx = self.line_of_sight(rx, edge_rx_coordinate)(distances)

        coordinate = np.argmin(np.abs(los_edge_tx-los_edge_rx))
        bullington_xaxis = distances[coordinate]
        bullington_yaxis = los_edge_rx[coordinate]

        self.elevation_plot.plot(bullington_xaxis, bullington_yaxis, "ko")
        self.elevation_plot.annotate("Bullington", xy=(bullington_xaxis, bullington_yaxis+5), color="k")
        return [bullington_xaxis, bullington_yaxis]

    def durkin(self, tx, rx, x ,y):
        edges = self.get_edge_points(tx, rx, x, y)
        complete_loss = self.loss(tx, rx, edges)
        return edges, complete_loss