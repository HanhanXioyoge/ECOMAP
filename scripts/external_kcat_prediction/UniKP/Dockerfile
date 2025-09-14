FROM python:3.7-slim-buster

COPY . .

RUN pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu113

RUN pip install -r requirements.txt