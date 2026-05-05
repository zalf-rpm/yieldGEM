soil_data = []
with open("Soil_coordinates_1km.csv", newline='') as _:
    for line in _.readlines()[1:]:
        es = line.split(",")
        soil_data.append({
            "id": int(es[0]),
            "filename": es[1],
            "lat": float(es[3]),
            "lon": float(es[2]),
            "field_id": float(es[4]),
            "sowing_day": float(es[5]),
            "nfertilizer": float(es[6]),
            "irrigation": float(es[7])
        })

print(soil_data)