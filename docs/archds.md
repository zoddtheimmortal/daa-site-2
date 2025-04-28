# Architecture and Dataset Info

## Architecture

All the datasets were run on the following machine:

```cpp
Hardware:
    Model Name: MacBook Pro
    Model Identifier: Mac15,6
    Model Number: MRX33HN/A
    Chip: Apple M3 Pro
    Total Number of Cores: 11 (5 performance and 6 efficiency)
    Memory: 18 GB
    System Firmware Version: 11881.81.4
    OS Loader Version: 11881.81.4

Software:
    System Version: macOS 15.3.2 (24D81)
    Kernel Version: Darwin 24.3.0
    Boot Volume: Macintosh HD
    Boot Mode: Normal
    Computer Name: rudy
    Secure Virtual Memory: Enabled
    System Integrity Protection: Enabled
    Time since boot: 21 days, 23 minutes
```

The machine was under medium load.

## Dataset Info

### Wikipedia vote network

The following dataset was extracted from [Wikipedia Vote Network](https://snap.stanford.edu/data/wiki-Vote.html).

-   This dataset includes complete administrator elections and voting history data extracted from Wikipedia’s edit history (up to **January 3, 2008**).
-   It consists of **2,794 elections** with a total of **103,663 votes**.
-   **7,066 users** participated, either by voting or being voted on

#### Dataset Statistics

-   **Nodes:** 7,115
-   **Edges:** 103,689
-   **Nodes in Largest WCC:** 7,066 _(0.993)_
-   **Edges in Largest WCC:** 103,663 _(1.000)_
-   **Nodes in Largest SCC:** 1,300 _(0.183)_
-   **Edges in Largest SCC:** 39,456 _(0.381)_
-   **Average Clustering Coefficient:** 0.1409
-   **Number of Triangles:** 608,389
-   **Fraction of Closed Triangles:** 0.04564
-   **Diameter (Longest Shortest Path):** 7
-   **90-Percentile Effective Diameter:** 3.8

### Enron Email Network

The following dataset was extracted from the [Enron Email Network](https://snap.stanford.edu/data/email-Enron.html).

-   This dataset covers email communication within a dataset of around **half a million emails**.
-   The data was made public by the **Federal Energy Regulatory Commission** during its investigation.
-   Nodes represent **email addresses**, and an undirected edge between nodes **i** and **j** indicates that **i** sent at least one email to **j**.
-   Non-Enron email addresses act as sinks and sources, meaning their communication is only observed with Enron addresses.

#### Dataset Statistics

-   **Nodes:** 36,692
-   **Edges:** 183,831
-   **Nodes in Largest WCC:** 33,696 _(0.918)_
-   **Edges in Largest WCC:** 180,811 _(0.984)_
-   **Nodes in Largest SCC:** 33,696 _(0.918)_
-   **Edges in Largest SCC:** 180,811 _(0.984)_
-   **Average Clustering Coefficient:** 0.4970
-   **Number of Triangles:** 727,044
-   **Fraction of Closed Triangles:** 0.03015
-   **Diameter (Longest Shortest Path):** 11
-   **90-Percentile Effective Diameter:** 4.8

### Skitter Topology Graph

The following dataset was extracted from the [CAIDA Skitter Internet Topology Dataset](http://www.caida.org/tools/measurement/skitter).

-   This dataset represents an **Internet topology graph** from traceroutes run daily in **2005**.
-   The data was collected from several scattered sources to millions of destinations.
-   It contains approximately **1.7 million nodes** and **11 million edges**.

#### Dataset Statistics

-   **Nodes:** 1,696,415
-   **Edges:** 11,095,298
-   **Nodes in Largest WCC:** 1,694,616 _(0.999)_
-   **Edges in Largest WCC:** 11,094,209 _(1.000)_
-   **Nodes in Largest SCC:** 1,694,616 _(0.999)_
-   **Edges in Largest SCC:** 11,094,209 _(1.000)_
-   **Average Clustering Coefficient:** 0.2581
-   **Number of Triangles:** 28,769,868
-   **Fraction of Closed Triangles:** 0.001802
-   **Diameter (Longest Shortest Path):** 25
-   **90-Percentile Effective Diameter:** 6
