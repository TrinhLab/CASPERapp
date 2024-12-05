FROM --platform=linux/amd64 python:3.11.9-slim

# Install system dependencies including X11, Qt dependencies, and required Linux libraries
RUN apt-get update && apt-get install -y \
    libgl1-mesa-glx \
    libx11-xcb1 \
    libxcb-icccm4 \
    libxcb-image0 \
    libxcb-keysyms1 \
    libxcb-randr0 \
    libxcb-render-util0 \
    libxcb-shape0 \
    libxcb-xfixes0 \
    libxcb-xinerama0 \
    libxkbcommon-x11-0 \
    xvfb \
    libegl1 \
    libopengl0 \
    libxcb-cursor0 \
    qt6-base-dev \
    glibc-source \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy the entire application
COPY . .

# Make sure the SeqFinder executable has correct permissions
RUN chmod +x /app/src/SeqFinder/Casper_Seq_Finder_Lin

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Set environment variables for Qt
ENV QT_QPA_PLATFORM=xcb
ENV XDG_RUNTIME_DIR=/tmp/runtime-root
ENV DISPLAY=:0

# Create runtime directory
RUN mkdir -p /tmp/runtime-root && chmod 0700 /tmp/runtime-root

# Command to run the application
CMD ["python3", "src/main.py"] 