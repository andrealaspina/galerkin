function [answer]=chatGPT(...
         question)
         % Ask chatGPT

% Import
import matlab.net.*;
import matlab.net.http.*;

% Define the API key
api_key='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';

% Define the API endpoint
api_endpoint='https://api.openai.com/v1/completions';

% Define the parameters for the API request
parameters=struct('prompt',question,'model','text-davinci-003','max_tokens',100);

% Define the headers for the API request
headers(1)=matlab.net.http.HeaderField('Content-Type','application/json');
headers(2)=matlab.net.http.HeaderField('Authorization',['Bearer ',api_key]);

% Define the request message
request=matlab.net.http.RequestMessage('post',headers,parameters);

% Send the request and store the response
response=send(request,URI(api_endpoint));

% Extract the answer (or the error)
if any(strcmp(fieldnames(response.Body.Data),'choices'))
  answer=response.Body.Data.choices(1).text;
else
  error(response.Body.Data.error.message);
end

% Remove initial newlines
while answer(1)==newline
  answer=answer(2:end);
end

% Add newline after each dot if answer is not a list
if not(contains(answer,'1. ') && contains(answer,'2. '))
  answer=replace(answer,'. ',['.',newline]);
end

% Diplay answer
fprintf('\n%s\n\n',answer);

end