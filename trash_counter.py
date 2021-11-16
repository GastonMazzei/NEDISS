
with open('tmp/trash.txt','r') as f:
        d = f.readlines()


results = {
    '[ANSWER_MESSAGES] recieved and answered':['We effectively captured',
        [0]],

    '[ANSWER_MESSAGES] saw but failed to catch':['We were faced with a probe that indicated ', [0]],

    '[ANSWER_MESSAGES] laps waiting for any':['Finished waiting one', [0]],

    '[ANSWER_MESSAGES] had to wait for sends to complete':['We are currently waiting for a send', [0]],

    "[ANSWER_MESSAGES] unfound nodes":[
        "[CRITICAL] THE NODE WAS NOT FOUND", [0]],


    '[PERFORM_REQUESTS] ix elements iterated':['SURVIVED 2',
        [0]],

    '[PERFORM_REQUESTS] elements ssent':['Ending blocking ssend',
        [0]],


    '[PERFORM_REQUESTS] thought it has been answered (12.500 & 0.000)':[
        'Recieved a response to our request',
        [0]],

    '[PERFORM_REQUESTS] obtained 12.5000':[
        'Recieved a response to our request: val 12',
        [0]],

    "[PERFORM_REQUESTS]  Total laps for Qpend.size() != target at perform_requests":[
        "About to compute 'QPend.size()", [0]],

    "[PERFORM_REQUESTS] sent message appears to have been arrived in perform_requests":["Case 2", [0]],

    "[PERFORM_REQUESTS] resends at perform_requests":["Case 1.2.3",[0]],

    "[PERFORM_REQUESTS] probe showed '1' so we started a reception in perform_requests":["Case 2.3.2", [0]],

    "[PERFORM_REQUESTS] recv has been cancelled in perform_requests":["Case 2.6.1",[0]],

    "[PERFORM_REQUESTS] recv existed and we computed its status":["Case 2.5",[0]],
    "[PERFORM_REQUESTS] recv existed, the status was not arrived so we flaged for recv restart":["Case 2.6.1",[0]],

    }

for _ in d:
    for k in results.keys():
        if results[k][0] in _:
            results[k][1][0]+=1;

for k in results.keys():
    print(f'{k} \t---\t---\t---\t {results[k][1]}\n')

