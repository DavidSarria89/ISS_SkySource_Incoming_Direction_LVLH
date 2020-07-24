* Python script to calculate the direction of a source in the sky (specified by Equatorial coordinates J2000) in the International Space Station (ISS) LVLH frame of reference. Also evaluates the minimum altitude reached by the light of sight (ISS<->Source). Uses Python Skyfield Library.
* Requires Python >=3.6
* Executes with `python Calculate_Sky_Source_incoming_dir_FULL.py` (or `python3 Calculate_Sky_Source_incoming_dir_FULL.py`) using the terminal.
* Edit `Calculate_Sky_Source_incoming_dir_FULL.py` (first lines) to change inputs.
* Inputs are a datetime (that is used to calculate ISS's position) and R.A. and declination (equatorial J2000) of the source in the sky.
* Maximum time Two-Line Element (TLE) for the ISS position calculation is 24/07/2020. Using larger dates is unreliable. TLE can be updated by changing the file `./dataFiles/ISS_orbital_info.txt`. Other satellites can also be added by using another TLE file.
