import customtkinter as ctk
from tkinter import *
from world import Map
from durkin import Durkin
from CTkMessagebox import CTkMessagebox


window = ctk.CTk()
window.title("Durkin Propagation Model")
window.geometry("750x750")
window.maxsize(width=750, height=750)
window.minsize(width=750, height=750)
ctk.set_default_color_theme("dark-blue")



def open_map():
    tx = [latitude_tx_entry, longitude_tx_entry]
    rx = [latitude_rx_entry, longitude_rx_entry]
    Map(tx=tx, rx=rx)

def durkin():
    try:
        tx = [float(longitude_tx_entry.get()), float(latitude_tx_entry.get())]
        rx = [float(longitude_rx_entry.get()), float(latitude_rx_entry.get())]
        f = float(frequency_entry.get())
        gt = float(txgain_entry.get())
        gr = float(rxgain_entry.get())
        htx = float(txheight_entry.get())
        hrx = float(rxheight_entry.get())
    except:
        CTkMessagebox(title="Error", message="Invalid parameters", icon="cancel")
    else:
        Durkin(tx,rx,f,gt,gr,htx,hrx)

tabview = ctk.CTkTabview(master=window, corner_radius=10, 
                         width=850, height=600)

insert_values = ctk.CTkLabel(master=window, text="Model parameters\t\t", font=("Arial Bold",25))
insert_values.grid(row=0, column=1, columnspan=2, padx=50, pady=20)

tx_label = ctk.CTkLabel(master=window, text="Tx:", font=("Arial",20))
tx_label.grid(row=1, column=0, padx=50, pady=20)

latitude_tx_entry = ctk.CTkEntry(master=window, width=250, height=40, font=("Arial",18), placeholder_text="latitude(degrees)")
latitude_tx_entry.grid(row=1, column=1, padx=10, pady=20)

longitude_tx_entry = ctk.CTkEntry(master=window, width=250, height=40, font=("Arial",18), placeholder_text="longitude(degrees)")
longitude_tx_entry.grid(row=1, column=2, padx=0, pady=20)


rx_label = ctk.CTkLabel(master=window, text="Rx:", font=("Arial",20))
rx_label.grid(row=2, column=0, padx=0, pady=10)

latitude_rx_entry = ctk.CTkEntry(master=window, width=250, height=40, font=("Arial",18), placeholder_text="latitude(degrees)")
latitude_rx_entry.grid(row=2, column=1, padx=10, pady=10)

longitude_rx_entry = ctk.CTkEntry(master=window, width=250, height=40, font=("Arial",18), placeholder_text="longitude(degrees)")
longitude_rx_entry.grid(row=2, column=2, padx=0, pady=10)


get_map = ctk.CTkButton(master=window, text="Get from Map", font=("Helvetica Bold",20), hover_color='#888f8e',
                             text_color="black", height=50, fg_color="#95948d", command=open_map)
get_map.grid(row=3, columnspan=2, column=1, padx=10, pady=10)


frequency_label = ctk.CTkLabel(master=window, text="Frequency:", font=("Arial",20))
frequency_label.grid(row=4, column=0, padx=20, pady=(50,20))

frequency_entry = ctk.CTkEntry(master=window, width=300, height=40, font=("Arial",18), placeholder_text="MHz")
frequency_entry.grid(row=4, column=1, padx=0, pady=(50,20))

txgain_label = ctk.CTkLabel(master=window, text="Tx Gain:", font=("Arial",20))
txgain_label.grid(row=5, column=0, padx=20, pady=20)

txgain_entry = ctk.CTkEntry(master=window, width=300, height=40, font=("Arial",18), placeholder_text="dimensionless")
txgain_entry.grid(row=5, column=1, padx=0)

rxgain_label = ctk.CTkLabel(master=window, text="Rx Gain:", font=("Arial",20))
rxgain_label.grid(row=6, column=0, padx=20, pady=20)

rxgain_entry = ctk.CTkEntry(master=window, width=300, height=40, font=("Arial",18), placeholder_text="dimensionless")
rxgain_entry.grid(row=6, column=1, padx=0, pady=20)

txheight_label = ctk.CTkLabel(master=window, text="Tx Height:", font=("Arial",20))
txheight_label.grid(row=7, column=0, padx=20, pady=20)

txheight_entry = ctk.CTkEntry(master=window, width=300, height=40, font=("Arial",18), placeholder_text="meters")
txheight_entry.grid(row=7, column=1, padx=0)

rxheight_label = ctk.CTkLabel(master=window, text="Rx Height:", font=("Arial",20))
rxheight_label.grid(row=8, column=0, padx=20, pady=20)

rxheight_entry = ctk.CTkEntry(master=window, width=300, height=40, font=("Arial",18), placeholder_text="meters")
rxheight_entry.grid(row=8, column=1, padx=0, pady=20)

calculate_button = ctk.CTkButton(master=window, text="Calculate", font=("Arial Bold",30), hover_color='#246629',
                             corner_radius=75, height=100, fg_color="#2A8C55", command=durkin)
calculate_button.grid(row=5, column=2, rowspan=3, padx=40)

window.mainloop()