//
// Created by m4zz31 on 31/10/21.
//

#ifndef GLOBAL_STANDARIZED_MESSAGES
#define GLOBAL_STANDARIZED_MESSAGES

const std::string msg_prev = "[info][informal-test] about to run ";
const std::string msg_post = "[info][informal-test] (apparent) success ";
const std::string min_threads = "[error] This program requires the max N of threads to be AT LEAST 2: master for message passing and one to be able to perform tasks. Environment variable is called OMP_THREAD_LIMIT, please set it to at least 2.";
const std::string queue_timeout = "[error] TIME_COUNTER_TOLERANCE has been exceded in Communication/CommunicationFunctions/GetAllMsgs. Queue was empty so Master timed out.";

#endif