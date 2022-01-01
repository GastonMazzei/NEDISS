<b>Index of the available tests</b>

* `graph-test-init.h`: <i> This test is for visually inspecting that the created graph has the correct number of nodes, edges and values. It uses a constant initialization.</i>

* `graph-test-singlestep-evolution.h`: <i>This test is for visually inspecting that the integration is yielding correct results with all the nodes and values initialized to the same value.</i>

* `long-singlestep-run.sh`: <i>This test initializes all nodes and edges with the same value, iterates the system several times, and finally shows the results. It aims to both, deepen in the previous test as to check that symmetric systems yield symmetric results, and also confirm that there is no memory leak nor segmentation fault that appears with a small probability.</i>

