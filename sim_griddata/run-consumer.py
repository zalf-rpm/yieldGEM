# #!/usr/bin/python
# # -*- coding: UTF-8

import csv
import json
import os
import sys
import zmq
import monica_io3


def run_consumer(server={"server": None, "port": None}):
    "collect data from workers"

    config = {
        "port": server["port"] if server["port"] else "7777",
        "server": server["server"] if server["server"] else "localhost",
        "path_to_out_file": os.path.join(os.path.dirname(__file__), "out.csv")
    }
    if len(sys.argv) > 1 and __name__ == "__main__":
        for arg in sys.argv[1:]:
            k, v = arg.split("=", maxsplit=1)
            if k in config:
                config[k] = v.lower() == "true" if v.lower() in ["true", "false"] else v
    print("config:", config)

    context = zmq.Context()
    socket = context.socket(zmq.PULL)

    socket.connect("tcp://" + config["server"] + ":" + config["port"])

    socket.RCVTIMEO = 6000  # stop consumer after 1 minutes
    leave = False

    def process_message(msg):

        if not hasattr(process_message, "received_env_count"):
            process_message.received_env_count = 0 # type: ignore
        if not hasattr(process_message, "header_written"):
            process_message.header_written = False # type: ignore

        leave = False

        if msg["type"] == "finish":
            print("c: received finish message")
            leave = True

        else:
            print("c: received work result ", process_message.received_env_count, " customId: ", # type: ignore
                  str(msg.get("customId", "")))

            process_message.received_env_count += 1 # type: ignore

            with open(config["path_to_out_file"], 'a', newline='') as _:

                writer = csv.writer(_, delimiter=",", lineterminator='\n')
                cid = msg["customId"]
                # print(msg)
                for data_ in msg.get("data", []):
                    results = data_.get("results", [])
                    print(results)

                    # Define output fields from your 1.json
                    # output_fields = [
                    #     "Month", "Year", "TraDef", "NDef", "Nstress",
                    #     "ActNup", "AbBiomNc", "LAI", "ETa/ETc", "PASW", "StomRes"
                    # ]
                    output_fields = ["Yield"]
                    # Write header row on first write
                    if not process_message.header_written: # type: ignore
                        header = ["ID"] + output_fields
                        writer.writerow(header)
                        process_message.header_written = True # type: ignore

                    # Extract data for each result
                    for result in results:
                        # Skip empty results
                        has_data = any(result.get(field) is not None and result.get(field) != ""
                                      for field in output_fields)
                        if not has_data:
                            continue

                        row = [cid["id"]]
                        for field in output_fields:
                            row.append(result.get(field, ""))
                        writer.writerow(row)


        return leave

    while not leave:
        try:
            msg = json.loads(socket.recv_string(encoding="latin-1"))
            leave = process_message(msg)
        except:
            print(sys.exc_info())
            leave = True
            # continue

    print("exiting run_consumer")


if __name__ == "__main__":
    run_consumer()