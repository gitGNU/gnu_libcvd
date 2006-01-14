/*                       
	This file is part of the CVD Library.

	Copyright (C) 2005 The Authors

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	This library is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
	Lesser General Public License for more details.

	You should have received a copy of the GNU Lesser General Public
	License along with this library; if not, write to the Free Software
	Foundation, Inc., 
    51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
// Paul Smith 1 March 2005
// Uses ffmpeg libraries to play most types of video file

#include <string>
#include <sstream>

#include <cvd/exceptions.h>
#include <cvd/videofilebuffer.h>

using namespace std;

namespace CVD
{

using namespace Exceptions::VideoFileBuffer;
	
//
// EXCEPTIONS
//

Exceptions::VideoFileBuffer::FileOpen::FileOpen(const std::string& name, const string& error)
{
	what = "RawVideoFileBuffer: Error opening file \"" + name + "\": " + error;
}

Exceptions::VideoFileBuffer::BadFrameAlloc::BadFrameAlloc()
{
	what = "RawVideoFileBuffer: Unable to allocate video frame.";
}

Exceptions::VideoFileBuffer::BadDecode::BadDecode(double t)
{
	ostringstream os;
	os << "RawVideoFileBuffer: Error decoding video frame at time " << t << ".";
	what = os.str();
}

Exceptions::VideoFileBuffer::EndOfFile::EndOfFile()
{
	what =  "RawVideoFileBuffer: Tried to read off the end of the file.";
}

Exceptions::VideoFileBuffer::BadSeek::BadSeek(double t)
{
	ostringstream ss;
	ss << "RawVideoFileBuffer: Seek to time " << t << "s failed.";
	what = ss.str();
}

namespace VFB
{
//
// CONSTRUCTOR
//
RawVideoFileBuffer::RawVideoFileBuffer(const std::string& file, bool rgbp) :
	end_of_buffer_behaviour(VideoBufferFlags::RepeatLastFrame),
	pFormatContext(0),
	pCodecContext(0),
	pFrame(0), 
	pFrameRGB(0),
	//buffer(0),
	frame_time(0.0),
	is_rgb(rgbp)
{
	try
	{
		// Register the formats and codecs
		av_register_all();
	
		// Now open the video file (and read the header, if present)
		if(av_open_input_file(&pFormatContext, file.c_str(), NULL, 0, NULL) != 0)
			throw FileOpen(file, "File could not be opened.");
		
		// Read the beginning of the file to get stream information (in case there is no header)
		if(av_find_stream_info(pFormatContext) < 0)
			throw FileOpen(file, "Stream information could not be read.");
		
		// Dump details of the video to standard error
		dump_format(pFormatContext, 0, file.c_str(), false);
		
		// We shall just use the first video stream
		video_stream = -1;
		for(int i=0; i < pFormatContext->nb_streams && video_stream == -1; i++)
		{
		    #if LIBAVFORMAT_BUILD >= 4629
			if(pFormatContext->streams[i]->codec->codec_type == CODEC_TYPE_VIDEO)
				video_stream = i; // Found one!
		    #else
			if(pFormatContext->streams[i]->codec.codec_type == CODEC_TYPE_VIDEO)
				video_stream = i; // Found one!
			#endif
		}
		if(video_stream == -1)
			throw FileOpen(file, "No video stream found.");
		
		// Get the codec context for this video stream
		#if LIBAVFORMAT_BUILD >= 4629
		pCodecContext = pFormatContext->streams[video_stream]->codec;
		#else
		pCodecContext = &pFormatContext->streams[video_stream]->codec;
		#endif
		
		// Find the decoder for the video stream
		AVCodec* pCodec = avcodec_find_decoder(pCodecContext->codec_id);
		if(pCodec == NULL)
		{
			pCodecContext = 0; // Since it's not been opened yet
			throw FileOpen(file, "No appropriate codec could be found.");
		}
		
		// Open codec
		if(avcodec_open(pCodecContext, pCodec) < 0)
		{
			pCodecContext = 0; // Since it's not been opened yet
			throw FileOpen(file, string(pCodec->name) + " codec could not be initialised.");
		}
		
		#if LIBAVCODEC_BUILD < 4754
		// Hack to fix wrong frame rates
		if(pCodecContext->frame_rate > 1000 && pCodecContext->frame_rate_base == 1)
			pCodecContext->frame_rate_base = 1000;
		#endif
		
		
		// Allocate video frame
		pFrame = avcodec_alloc_frame();
		if(pFrame == NULL)
			throw BadFrameAlloc();
		
		// And a frame to hold the RGB version
		pFrameRGB = avcodec_alloc_frame();
		if(pFrameRGB == NULL)
			throw BadFrameAlloc();
		
		// How big is the buffer?
		//long num_bytes = avpicture_get_size(PIX_FMT_RGB24, pCodecContext->width, pCodecContext->height);
		my_size = ImageRef(pCodecContext->width, pCodecContext->height);

		// And allocate a contiguous buffer
		//buffer = new CVD::Rgb<CVD::byte>[my_size.x * my_size.y];

		// Assign this buffer to image planes in pFrameRGB
		//avpicture_fill((AVPicture *)pFrameRGB, reinterpret_cast<uint8_t*>(buffer), PIX_FMT_RGB24, pCodecContext->width, pCodecContext->height);
	
		// Now read the first frame
		if(!read_next_frame())
			throw EndOfFile();
		
		start_time = 0;
		frame_ready = true;
	}
	catch(CVD::Exceptions::All)
	{
		// Tidy things up on the heap if we failed part-way through constructing
		if(pFormatContext != 0)
		    av_close_input_file(pFormatContext);
		
		if(pCodecContext != 0)
			avcodec_close(pCodecContext);
		
		if(pFrame != 0)
			av_free(pFrame);
		
		if(pFrameRGB != 0)
			av_free(pFrameRGB);
		
		//if(buffer != 0)
		//	delete[] buffer;
		
		// Now re-throw
		throw;
	}
}

//
// DESTRUCTOR
//
RawVideoFileBuffer::~RawVideoFileBuffer()
{
    //delete [] buffer;
    av_free(pFrameRGB);
    av_free(pFrame);
    avcodec_close(pCodecContext);
    av_close_input_file(pFormatContext);
}


//
// READ NEXT FRAME
//
bool RawVideoFileBuffer::read_next_frame()
{
	//Make next_frame point to a new block of data, getting the sizes correct.
	if(is_rgb)
	{
		Image<Rgb<byte> > tmp(my_size);
		next_frame = (reinterpret_cast<Image<byte>&>(tmp));
	}
	else
	{
		Image<byte> tmp(my_size);
		next_frame = tmp;
	}

	//Assign this new memory block 
	avpicture_fill((AVPicture *)pFrameRGB, reinterpret_cast<uint8_t*>(next_frame.data()), is_rgb?PIX_FMT_RGB24:PIX_FMT_GRAY8, pCodecContext->width, pCodecContext->height);


    AVPacket packet;
	packet.stream_index = -1;
	
	// How many frames do we read looking for our video stream?
	// If we assume our streams are interlaced, and some might be interlaced
	// 2:1, this should probably do
	const int max_loop = MAX_STREAMS * 2; 
	
	int i;
	for(i = 0; packet.stream_index != video_stream && i < max_loop; i++)
	{
		
    	if(av_read_frame(pFormatContext, &packet) < 0)
			return false;

		if(packet.stream_index == video_stream)
		{
			// Ask this packet what time it is
			if(packet.pts >= 0)
				frame_time = packet.pts / static_cast<double>(AV_TIME_BASE);
			else // sometimes this is reported incorrectly, so guess
			{
				frame_time = frame_time + 1.0 / frames_per_second();
			}
			
			// Decode video frame
			int got_picture;
			if(avcodec_decode_video(pCodecContext, pFrame, &got_picture, 
				packet.data, packet.size) == -1)
			{
				throw BadDecode(frame_time);
			}
	
			// Did we get a video frame?
			if(got_picture)
			{
				// Convert the image from its native format to RGB
				img_convert((AVPicture *)pFrameRGB, is_rgb?PIX_FMT_RGB24:PIX_FMT_GRAY8, 
					(AVPicture*)pFrame, pCodecContext->pix_fmt, 
					pCodecContext->width, pCodecContext->height);
				
			}
		}

		// Free the packet that was allocated by av_read_frame
		av_free_packet(&packet);
	}
	
	// Did we not find one?
	if(i == max_loop)
		return false;
	
	return true;
}


//
// GET FRAME
//
VideoFileFrame<byte>* RawVideoFileBuffer::get_frame()
{

	if(!frame_pending())
		throw EndOfFile();

// 	Don't use - pCC->frame_number doesn't reset after a seek!
//  Instead, we ask the packet its time when we decode it
//	double time = start_time + pCodecContext->frame_number * pCodecContext->frame_rate_base / static_cast<double>(pCodecContext->frame_rate);
	VideoFileFrame<byte>* vf = new VideoFileFrame<byte>(frame_time, next_frame);

	if(!read_next_frame())
	{
		switch(end_of_buffer_behaviour)
		{
			case VideoBufferFlags::RepeatLastFrame:
				// next_frame is empty because there isn't one, so 
				// I'll copy the one that I'm about to return so that
				// I can return it next time as well
				if(is_rgb)
				{
					Image<Rgb<byte> > tmp = reinterpret_cast<Image<Rgb<byte> >&>(next_frame);
					tmp.copy_from(reinterpret_cast<VideoFileFrame<Rgb<byte> >&>(*vf));
					next_frame = (reinterpret_cast<Image<byte>&>(tmp));
				}
				else
				{
					next_frame.copy_from(*vf);
				}
				break;
			
			case VideoBufferFlags::UnsetPending:
				frame_ready = false;
			   break;
			
			case VideoBufferFlags::Loop:
				seek_to(start_time);
				break;
		}
	}
	
	return vf;	
}

//
// PUT FRAME
//
void RawVideoFileBuffer::put_frame(VideoFrame<byte>* f)
{
	VideoFileFrame<byte>* vff  = dynamic_cast<VideoFileFrame<byte> *>(f);

	if(!vff)
		throw Exceptions::VideoBuffer::BadPutFrame();
	else
		delete vff;
}

//
// SEEK TO
//
void RawVideoFileBuffer::seek_to(double t)
{	
	#if LIBAVFORMAT_BUILD >= 4623
	if(av_seek_frame(pFormatContext, -1, static_cast<int64_t>(t*AV_TIME_BASE+0.5), AVSEEK_FLAG_ANY) < 0)
	#else
	if(av_seek_frame(pFormatContext, -1, static_cast<int64_t>(t*AV_TIME_BASE+0.5)) < 0)
	#endif
	{
		cerr << "av_seek_frame not supported by this codec: performing (slow) manual seek" << endl;
		
		// Seeking is not properly sorted with some codecs
		// Fudge it by closing the file and starting again, stepping through the frames
		string file = pFormatContext->filename;
		av_close_input_file(pFormatContext);
		avcodec_close(pCodecContext);
		
		// Now open the video file (and read the header, if present)
		if(av_open_input_file(&pFormatContext, file.c_str(), NULL, 0, NULL) != 0)
			throw FileOpen(file, "File could not be opened.");
		
		// Read the beginning of the file to get stream information (in case there is no header)
		if(av_find_stream_info(pFormatContext) < 0)
			throw FileOpen(file, "Stream information could not be read.");
		
		// No need to find the stream--we know which one it is (in video_stream)
		
		// Get the codec context for this video stream
		#if LIBAVFORMAT_BUILD >= 4629
		pCodecContext = pFormatContext->streams[video_stream]->codec;
		#else
		pCodecContext = &pFormatContext->streams[video_stream]->codec;
		#endif
		
		// Find the decoder for the video stream
		AVCodec* pCodec = avcodec_find_decoder(pCodecContext->codec_id);
		if(pCodec == NULL)
		{
			pCodecContext = 0; // Since it's not been opened yet
			throw FileOpen(file, "No appropriate codec could be found.");
		}
		
		// Open codec
		if(avcodec_open(pCodecContext, pCodec) < 0)
		{
			pCodecContext = 0; // Since it's not been opened yet
			throw FileOpen(file, string(pCodec->name) + " codec could not be initialised.");
		}
		
		start_time = 0;
		frame_ready = true;
		
		// REOPENED FILE OK
		// Now read frames until we get to the time we want
		
		int frames = static_cast<int>((t * frames_per_second() + 0.5));
		for(int i = 0; i < frames; i++)
		{
			read_next_frame();
		}
	}
	
	if(!read_next_frame())
			throw BadSeek(t);
}

}
} // namespace CVD