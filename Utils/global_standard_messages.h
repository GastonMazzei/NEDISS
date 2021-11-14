//
// Created by m4zz31 on 31/10/21.
//

#ifndef GLOBAL_STANDARIZED_MESSAGES
#define GLOBAL_STANDARIZED_MESSAGES
#include <string>

const std::string separator = "-------------------------------------";
const std::string msg_prev = "[info][informal-test] about to run ";
const std::string msg_post = "[info][informal-test] (apparent) success ";
const std::string min_threads = "[error] This program requires the max N of threads to be AT LEAST 3: 1 for master (message passing), 1 for the master assistant and 1 for the worker. Environment variable is called OMP_THREAD_LIMIT, please set it to at least 3.";
const std::string nested_threads = "[error] This program requires the OMP_NESTED be set to true. Master needs an assistant!";
const std::string queue_timeout = "[error] TIME_COUNTER_TOLERANCE has been exceded in Communication/CommunicationFunctions/GetAllMsgs. Queue was empty so Master timed out.";
const std::string init_warning = "[warning] We currently rely heavily on the long -> int conversion which during compilation raises warnings of being implementation-defined. Two important lines are:    \n(1) 'rank = in_degree(*v, g) + out_degree(*v, g);' and\n(2) 'NLocals = neighbors.second - neighbors.first;'";
const std::string LOGO = "************************************************\n"
                         "     .-\"\"L_ / \\  /|/  __//  _ \\/ \\/ ___\\/ ___\\  \n"
                         ";`, /   ( o\\| |\\ |||  \\  | | \\|| ||    \\|    \\  \n"
                         "\\  ;    `, /| | \\|||  /_ | |_/|| |\\___ |\\___ |  \n"
                         ";_/\"`.__.-\" \\_/  \\|\\____\\\\____/\\_/\\____/\\____/  \n"
                         "************************************************\n"
                         "  Network Diffusion & Synchronization Simulator";
const std::string init_msg = "\n\n"+LOGO+"\n\n\n";
const std::string min_subthread_msg = "[error] Minimum Required Subthreads are at least 1.";
const std::string min_batch_msg = "[error] a batch size of at least 1 is required in the message handling interface.";
const std::string architecture_and_maxtag = "[error] ";

#endif