# BlockDrift: Exact Top-k Changes Between Historical Windows from Paired Prefix Bounds

BlockDrift ranks keys by the magnitude of their signed count change between two historical windows. Two endpoint pairings bound each block’s change. Exact counts settle the surviving ranks and ties.

## Run

64-bit Windows, Python 3.10+ and Visual Studio C++ Build Tools.

~~~powershell
.\build.cmd
python run_example.py
~~~

[Data](data) includes the window example and paired service and construction costs for Chicago 311, New York taxis, Divvy and Citi Bike.
