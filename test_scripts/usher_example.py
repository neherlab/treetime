import gzip
import json
from google.protobuf import text_format
from google.protobuf.json_format import MessageToDict

from usher import parsimony_pb2


def read_parsimony_data_pb_gz(file_path):
  parsimony_data = parsimony_pb2.data()
  with gzip.open(file_path, 'rb') as f:
    parsimony_data.ParseFromString(f.read())
  return parsimony_data


def write_parsimony_data_pb_gz(file_path, data):
  parsimony_data = parsimony_pb2.data()
  parsimony_data.CopyFrom(data)
  with gzip.open(file_path, 'wb') as f:
    f.write(parsimony_data.SerializeToString())


def read_parsimony_pb_text(file_path):
  parsimony_data = parsimony_pb2.data()
  with open(file_path, 'r') as f:
    text_format.Merge(f.read(), parsimony_data)
  return parsimony_data


def write_parsimony_pb_text(file_path, data):
  parsimony_data = parsimony_pb2.data()
  parsimony_data.CopyFrom(data)
  with open(file_path, 'w') as f:
    f.write(text_format.MessageToString(parsimony_data))


def main():
  # Read `parsimony_pb2.data` message from gzipped binary file
  data = read_parsimony_data_pb_gz("data/usher/latest/public-latest.all.masked.msa.pb.gz")

  # Convert protobuf object to dict
  data_dict = MessageToDict(data)

  # Convert dict to json string
  data_json = json.dumps(data_dict, indent=2)
  print(data_json)

  # Write original `parsimony_pb2.data` message into "text format" file
  write_parsimony_pb_text("data/usher/latest/public-latest.all.masked.msa.pbtxt", data)


if __name__ == '__main__':
  main()
