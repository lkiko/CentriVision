# -*- coding: utf-8 -*-
"""
MP4 转压缩 GIF（兼容 moviepy 2.x + Pillow 10+ + NumPy 2.x）
"""
# # # conda activate py312

from moviepy.video.io.VideoFileClip import VideoFileClip
from PIL import Image

def mp4_to_gif(input_path, output_path, width=320, fps=10, duration=None):
    """
    将 MP4 转为压缩 GIF
    :param input_path: 输入 MP4 文件路径
    :param output_path: 输出 GIF 文件路径
    :param width: GIF 宽度，保持纵横比
    :param fps: 帧率，越低越小
    :param duration: 可选，裁剪视频前几秒 (start, end) 或 None 全部
    """
    # 兼容旧 MoviePy 代码使用 ANTIALIAS
    if not hasattr(Image, "ANTIALIAS"):
        Image.ANTIALIAS = Image.Resampling.LANCZOS

    # 加载视频
    clip = VideoFileClip(input_path)

    # 裁剪时长
    if duration:
        if isinstance(duration, tuple) and len(duration) == 2:
            clip = clip.subclipped(duration[0], duration[1])
        else:
            raise ValueError("duration 参数应为 (start, end) 或 None")

    # 缩小尺寸
    clip_resized = clip.resized(width=width)

    # 输出 GIF（MoviePy 2.x 不再支持 program 参数）
    clip_resized.write_gif(output_path, fps=fps)
    clip.close()
    print(f"GIF 已生成: {output_path}")


# -----------------------------
# 使用示例
# -----------------------------
if __name__ == "__main__":
    input_mp4 = "install.mp4"
    output_gif = "install.gif"
    # 转换前 5 秒，宽度 320，帧率 10
    mp4_to_gif(input_mp4, output_gif, width=1000, fps=10, duration=(0,5))






