import customtkinter
from tkinter import *
from tkintermapview import TkinterMapView
from CTkMessagebox import CTkMessagebox

customtkinter.set_default_color_theme("dark-blue")


class Map(customtkinter.CTkToplevel):
    
    APP_NAME = "Position Selector"
    WIDTH = 850
    HEIGHT = 600

    def __init__(self, tx ,rx, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.tx = tx
        self.rx = rx
        self.title(Map.APP_NAME)
        self.geometry(str(Map.WIDTH) + "x" + str(Map.HEIGHT))
        self.minsize(Map.WIDTH, Map.HEIGHT)

        self.marker_list = dict.fromkeys(["tx", "rx", "tx_marker", "rx_marker"])

        # ============ create two CTkFrames ============

        self.grid_columnconfigure(0, weight=0)
        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=1)

        self.frame_left = customtkinter.CTkFrame(master=self, width=150, corner_radius=0, fg_color=None)
        self.frame_left.grid(row=0, column=0, padx=0, pady=0, sticky="nsew")

        self.frame_right = customtkinter.CTkFrame(master=self, corner_radius=0)
        self.frame_right.grid(row=0, column=1, rowspan=1, pady=0, padx=0, sticky="nsew")

        # ============ frame_left ============

        self.frame_left.grid_rowconfigure(2, weight=1)

        self.confirm = customtkinter.CTkButton(master=self.frame_left,
                                                text="Confirm Positions",
                                                command=self.confirm_markers)
        self.confirm.grid(pady=(100, 0), padx=(20, 20), row=0, column=0)

        self.clear = customtkinter.CTkButton(master=self.frame_left,
                                                text="Clear",
                                                command=self.clear_marker_event)
        self.clear.grid(pady=(20, 0), padx=(20, 20), row=1, column=0)

        self.map_label = customtkinter.CTkLabel(self.frame_left, text="Map:", anchor="w")
        self.map_label.grid(row=3, column=0, padx=(20, 20), pady=(20, 0))
        self.map_option_menu = customtkinter.CTkOptionMenu(self.frame_left, values=["OpenStreetMap", "Google", "Google Satellite"],
                                                                       command=self.change_map)
        self.map_option_menu.grid(row=4, column=0, padx=(20, 20), pady=(10, 100))

        # ============ frame_right ============

        self.frame_right.grid_rowconfigure(1, weight=1)
        self.frame_right.grid_rowconfigure(0, weight=0)
        self.frame_right.grid_columnconfigure(0, weight=1)
        self.frame_right.grid_columnconfigure(1, weight=0)
        self.frame_right.grid_columnconfigure(2, weight=1)

        self.world = TkinterMapView(self.frame_right, corner_radius=0)
        self.world.grid(row=1, rowspan=1, column=0, columnspan=3, sticky="nswe", padx=(0, 0), pady=(0, 0))

        # Set default values
        self.world.set_zoom(0)
        self.map_option_menu.set("OpenStreetMap")
        self.world.add_right_click_menu_command(label="Add Tx",
                                command=self.add_tx,
                                pass_coords=True)
        self.world.add_right_click_menu_command(label="Add Rx",
                                command=self.add_rx,
                                pass_coords=True)

    def confirm_markers(self):
        if self.marker_list["tx"] and self.marker_list["rx"]:
            self.tx[0].delete(0,END)
            self.tx[1].delete(0,END)
            self.rx[0].delete(0,END)
            self.rx[1].delete(0,END)
            self.tx[0].insert(END, self.marker_list["tx"][0])
            self.tx[1].insert(END, self.marker_list["tx"][1])
            self.rx[0].insert(END, self.marker_list["rx"][0])
            self.rx[1].insert(END, self.marker_list["rx"][1])

            self.destroy()
        else:
            CTkMessagebox(title="Error", message="Insert Tx and Rx", icon="cancel")

    def clear_marker_event(self):
        if self.marker_list["tx_marker"]:
            self.marker_list["tx_marker"].delete()
        if self.marker_list["rx_marker"]:
            self.marker_list["rx_marker"].delete()

    def change_map(self, new_map: str):
        if new_map == "OpenStreetMap":
            self.world.set_tile_server("https://a.tile.openstreetmap.org/{z}/{x}/{y}.png")
        elif new_map == "Google":
            self.world.set_tile_server("https://mt0.google.com/vt/lyrs=m&hl=en&x={x}&y={y}&z={z}&s=Ga", max_zoom=22)
        elif new_map == "Google Satellite":
            self.world.set_tile_server("https://mt0.google.com/vt/lyrs=s&hl=en&x={x}&y={y}&z={z}&s=Ga", max_zoom=22)

    def add_tx(self, coords: float):
        print("Tx marker:", coords)
        if self.marker_list["tx_marker"]:
            self.marker_list["tx_marker"].delete()
        marker = self.world.set_marker(coords[0], coords[1], 
            text="Transmitter\n({:.2f}, {:.2f})".format(coords[0], coords[1]),
            marker_color_circle="#0000FF", marker_color_outside="#428df5", 
            font=("Arial Bold",12), text_color="black")
        self.marker_list["tx"] = coords
        self.marker_list["tx_marker"] = marker

    def add_rx(self, coords: float):
        print("Rx marker:", coords)
        if self.marker_list["rx_marker"]:
            self.marker_list["rx_marker"].delete()
        marker =self.world.set_marker(coords[0], coords[1], 
            text="Receiver\n({:.2f}, {:.2f})".format(coords[0], coords[1]),
            marker_color_circle="#32a885", marker_color_outside="#5fe18d",
            font=("Arial Bold",12), text_color="black")
        self.marker_list["rx"] = coords
        self.marker_list["rx_marker"] = marker