from PIL import Image
import os

# 連番の画像ファイルが保存されているディレクトリとファイル名のパターンを指定
image_folder = '/Users/takinosenaoto/Downloads/my_research/CS用/動画'  # 画像ファイルのフォルダ
output_video = '/Users/takinosenaoto/Downloads/my_research/CS用/動画/output_video2.gif'  # 出力ファイル名

# 画像ファイルの拡張子と連番の開始番号を指定
file_extension = 'png'
start_frame = 0

# 画像ファイルのリストを作成
images = []
for i in range(start_frame, 10000):  # 必要に応じて上限値を調整
    filename = f"{image_folder}/output2.{str(i).zfill(4)}.{file_extension}"
    if os.path.isfile(filename):
        images.append(Image.open(filename))

# GIFとして保存
if images:
    images[0].save(output_video, save_all=True, append_images=images[1:], duration=40, loop=0)

